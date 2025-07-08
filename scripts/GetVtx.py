#!/usr/bin/python3
#
# This is a program for getting versions of 'bossver' which can be
# used to query information about the vertex from the database.
#
# GetVtx(run) - print a list of available BossVer for a given run
#
# GetSftVer(runbeg,runend) - print SftVer(BossVer) for given range of
#                            runs and given BossRelease

import MySQLdb

class RunInfo:
    # ----------------------------------------------------------------
    def __init__(self):
    # ----------------------------------------------------------------
        self.verbose = False    # default: silent
        self.bossver = "6.6.4"
        self.host = "bes3db2.ihep.ac.cn"
        self.offlinedb = MySQLdb.connect (host =  self.host,
                                          user = "guest",
                                          passwd = "guestpass",
                                          db = "offlinedb")

    # print Boss version(s) for a given run:
    # ----------------------------------------------------------------
    def GetVtx(self,run):
    # ----------------------------------------------------------------
        cur_offline = self.offlinedb.cursor()
        sql = "SELECT Vx, Vy, Vz, SftVer"\
              " FROM BeamPar"\
              " WHERE `RunNo`='{RunNo}'"
        sql = sql.format(RunNo=run)
        if self.verbose:
            print(sql)
        cur_offline.execute (sql)
        if self.verbose:
            print(" count= ",cur_offline.rowcount)
        print("run#",run,":")
        if cur_offline.rowcount == 0:
            print("ERROR: can not find records")
        else:
            so="     BossVer= {:15s} Vtx=[{:.4f},{:.4f},{:.4f}]"
            for row in cur_offline.fetchall():
                (Vx,Vy,Vz,SftVer)=row
                print(so.format(SftVer,Vx,Vy,Vz))
        cur_offline.close()

    # print SftVer for given range of runs and given BossRelease
    # ----------------------------------------------------------------
    def GetSftVer(self,runbeg,runend):
    # ----------------------------------------------------------------
        cur_offline = self.offlinedb.cursor()
        sql = "SELECT SftVer, ParVer"\
              " FROM CalVtxLumVer"\
              " WHERE BossRelease='{BossRelease}' and"\
              " RunFrom <= '{RunFrom}' and RunTo >= '{RunTo}'"\
              " and DataType='LumVtx'"
        sql = sql.format(BossRelease=self.bossver,
                         RunFrom=runbeg,RunTo=runend)
        if self.verbose:
            print(sql)
        cur_offline.execute (sql)
        print("runs[",runbeg,runend,"] -> ", end=' ')
        if cur_offline.rowcount == 0:
            print("ERROR: can not find records")
        else:
            (SftVer,ParVer)=cur_offline.fetchone()
            print("BossVer=",SftVer," ParVer=",ParVer)
        cur_offline.close()

# Examples to use ----------------------------------------------------
r = RunInfo()
#  r.verbose = True

# phasejpsiscan
#  r.GetVtx(28312)
#  r.bossver = "6.6.4"
#  r.GetSftVer(28312, 28346) # 3050

#  r.GetVtx(28241)
#  r.GetVtx(27147)
#  r.bossver = "6.6.3"
#  r.GetSftVer(28241, 28266) # 3080 old
#  r.GetSftVer(27147,27233)  # 3082 old 3080 - I

# J/Psi 2009
#  r.GetVtx(9947)
#  r.GetVtx(10878)
#  r.GetVtx(10001)
#  r.bossver = "6.6.3"
#  r.GetSftVer(9947, 10878)

# R-scan 2015
#  r.GetVtx(39355)
#  r.bossver = "6.6.5.p01"
#  r.bossver = "7.0.3"
#  r.GetSftVer(39355, 39618) # 3080 new

# J/Psi scan 2018
#  r.GetVtx(55060)
#  r.bossver = "7.0.4"
#  r.GetSftVer(55060, 55065) # J1

# 3080 data 2019
#  r.GetVtx(59016)
#  r.bossver = "7.0.4"
#  r.GetSftVer(59016, 59141)

# Psi(2S) 2009 : run8093-run9025
r.GetVtx(8093)
#  r.bossver = "7.0.9"
#  r.GetSftVer(8093,9025)

# Psi(2S) 2012: run25338-run27090
#  r.GetVtx(25338)
#  r.GetSftVer(25338,27090)

# Psi(2S) 2021: run66257-run69292
#  r.GetVtx(66257)
#  r.bossver = "7.0.9"
#  r.GetSftVer(66257,69292)

# 3650 data: run9613-9779 , run33725-33772, run69612-run70132
#  r.GetVtx(9613)
#  r.GetVtx(33725)
#  r.GetVtx(69612)

