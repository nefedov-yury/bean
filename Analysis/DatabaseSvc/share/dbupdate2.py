#!/usr/bin/env python2
import sys
#import hashlib
#import md5
import os, os.path
import tempfile
import urllib2
import time

#  HOST = "http://bes3.jinr.ru/db_ana/"
HOST = "http://docbes3.ihep.ac.cn/db_ana/"
BLOCKSIZE = 1024*1024
gzsuffix = ".gz"
tmpsuffix=".tmp"

def hexify(s):
    return ("%02x"*len(s)) % tuple(map(ord, s))

def md5sum(path):
        sumstring = ''
        try:
                f = open(path,'rb')
                try:
                        import hashlib
                        sum = hashlib.md5()
                except ImportError:
                        # for Python << 2.5
                        import md5
                        sum = md5.new()

                while 1:
                        block = f.read(BLOCKSIZE)
                        if not block:
                                break
                        sum.update(block)
                sumstring = hexify(sum.digest())
                f.close()
        except:
                print "Error computing md5checksum of ["+path+"]."

        return sumstring


def CheckDirectory(destdir):
    if(not os.path.exists(destdir)):
        os.mkdir(destdir)

    if (os.path.exists(destdir) and not os.path.isdir(destdir)):
        print destdir+" exists, and it is not a directory!"
        return False

    if (not os.access(destdir,os.F_OK and os.W_OK and os.X_OK)):
        print destdir + " is not writeable"

    os.chdir(destdir)
    return True

def IsLocked():
    try:
     urllib2.urlopen(HOST+"/LOCK")
     return True
    except urllib2.URLError,e:
     msg = "Can't access server "+HOST
     if hasattr(e, 'reason'):
     #some common e.reason values:
     #       (-2, 'Name or service not known')
     #       (113, 'No route to host')
        try: msg += ": [%s]" % e.reason[1]
        except IndexError: pass
     if hasattr(e, 'code') and hasattr(e, 'msg'):
     #in this case e is HTTPError, a subclass of URLError)
     #some common *e.code, e.msg) values:
     #       (404, 'Not Found')
     #       (403, 'Forbidden')
     #       (500, 'Internal Server Error')
        if(e.code==404): return False;
        msg += ": HTTP error [%s, %s]" % (e.code, e.msg)
        msg += "."
     print msg
     exit(-1)

def GetTimestamp():
    Files={}
    try:
     opener = urllib2.build_opener()
     sock = opener.open(HOST+"/db.timestamp")
     text = sock.read()
     sock.close()
     lines=text.splitlines()
    except urllib2.URLError,e:
     msg = "Can't access server "+HOST
     if hasattr(e, 'reason'):
     #some common e.reason values:
     #       (-2, 'Name or service not known')
     #       (113, 'No route to host')
        try: msg += ": [%s]" % e.reason[1]
        except IndexError: pass
     if hasattr(e, 'code') and hasattr(e, 'msg'):
     #in this case e is HTTPError, a subclass of URLError)
     #some common *e.code, e.msg) values:
     #       (404, 'Not Found')
     #       (403, 'Forbidden')
     #       (500, 'Internal Server Error')
        msg += ": HTTP error [%s, %s]" % (e.code, e.msg)
     print msg
     exit(-1)
    except:
     print "Could not retrieve timestamp"
     exit(-1)
    try:
      timestamp=int(lines[0])
    except:
      print lines
      exit(-1)
    for line in lines[1:]:
      temp1,temp2 = line.split()
      Files[temp2]=temp1
    return timestamp,Files

def SelectDbForUpdate(timestamp,Files):
    LoadList = []

    if(os.path.exists("db.timestamp")):
      ts = open("db.timestamp")
      lines = ts.readlines();
      ts.close()

      if(not lines[0].strip().isdigit()):
        print "File db.timestamp is broken. Re-installing everything from scratch."
        try: os.unlink("db.timestamp")
        except: pass
        for fname in Files.keys():
         LoadList.append(fname)
        return LoadList

      cur_timestamp = int(lines[0])

      if( cur_timestamp == timestamp ):
       print "Warning: current database is up-to-date. Do nothing."
       return LoadList
      elif( cur_timestamp >= timestamp ):
       print "Warning! Please check: current database is NEWER than the remote replica. Do nothing."
       return LoadList

      for fname in Files.keys():
        if(not os.path.exists(fname)):
          print "Warning: " + fname + " is lost"
          LoadList.append(fname)
        else:
          chksum = md5sum(fname)
          if(not Files[fname] == chksum):
             print "Checksum of "+fname+" is different"
             LoadList.append(fname)
    else:
      print "No db.timestamp is present. Considering fresh installation."
      for fname in Files.keys():
         LoadList.append(fname)
    return LoadList


def DownloadDb(dbname):
    tmpdbname = dbname+tmpsuffix+gzsuffix
    remotedbname = dbname+gzsuffix
    if( os.path.exists(tmpdbname)):
      os.unlink(tmpdbname)

    try:
       # check lock
       for i in range(10):
        if(IsLocked()):
          var=""
          if(i>0): var="\r"
          sys.stdout.write(var+"Lock file found. It looks like server replica is being updated. Sleeping for "+str(i+1)+" min")
          sys.stdout.flush()
          time.sleep(1)
        elif (i>0):
          print "\nThe lock was removed"
          break
        else:
          break
       if(i==9):
          print "\nTimeout: 10 min waiting to unlock the server. Please try again later, or complain."
          exit(-1)

       ret = os.system("wget "+HOST+"/"+remotedbname+" -O "+tmpdbname)
    except os.error,e:
       print "Error %d: %s" % (e.args[0], e.args[1])
       try:os.unlink(tmpdbname)
       except: pass
       return ""
    except:
       print "Unable download file "+remotedbname
       try: os.unlink(tmpdbname)
       except: pass
       exit(-1)

    if(ret != 0 or not os.path.exists(tmpdbname) or os.path.getsize(tmpdbname)==0):
       print "Unable download database "+tmpdbname
       if(os.path.exists(tmpdbname)):
         os.unlink(tmpdbname)
       return ""

    ret = os.system("gunzip -f "+tmpdbname)
    if(not ret == 0):
       print "Unable to uncompress file "+tmpdbname
       return ""

    return dbname+tmpsuffix

def main(argv):
    argc=len(argv)
    destdir=""
    if(argc==2):
       destdir=argv[1]
    elif(argc>2):
       print "Usage: dbupdate.py [LOCAL_DB_DIR]"
       exit(0)

    if(destdir==""):
      var = os.getenv("SITEROOT")
      if(var <> None):
        destdir=os.path.join(var,"database")

    if(destdir==""):
      print "Destination directory neither given explicitly, nor can be derived from $BesArea. Use local directory."
      destdir="./"

    print "Install/update databases in the directory "+destdir

# Check wget is in the PATH
    if(os.system("wget > /dev/null 2>&1")==32512):
      print "Unable to find wget. Please install it first"
      exit(-1);

# Change directory to destdir and check if it is writable
    if(not CheckDirectory(destdir)):
      exit(-1)

# Get timestamp and checksums
    timestamp,Files = GetTimestamp()

# Check if db files are already loaded, and if they should be updated
    ListToLoad = SelectDbForUpdate(timestamp,Files)

    if(len(ListToLoad)>0):
       print "Going to download files: ", ListToLoad

       # Download databases to temporary files
       ListDone = []
       for db in ListToLoad:
         fname=DownloadDb(db)
         if(not fname == ""):
           ListDone.append(fname)

       # Check size and md5sum and rename temporary db to production db
       for tmpdbname in ListDone:
          chksum = md5sum(tmpdbname)
          dbname = tmpdbname.replace(tmpsuffix,"")
          if(chksum == Files[dbname]):
             os.rename(tmpdbname,dbname)
             print "Update database "+dbname
          else:
             print "Error: MD5 checksum is wrong for file "+tmpdbname

       # Update timestamp
       f=open("db.timestamp","w")
       f.truncate()
       f.write(str(timestamp)+"\n")
       for db in Files.keys():
         f.write(Files[db]+" "+db+"\n")
       f.close()

# Exit
#    exit(0)

if __name__ == '__main__':
    main(sys.argv)

