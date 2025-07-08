// Functions to copy files and directories. This code is highly
// operating system dependent and has been tested on Linux and MacOS.
// TODO: use C++17 filesystem library
//
// The BEAN uses these functions to copy the sqlite database
// to shared memory or /tmp.

#include <cstdio>
#include <cstdlib>
#if defined (__APPLE__)
#include <unistd.h>   // for mkdtemp()
#endif
#include <cerrno>

// opendir,mkdir,closedir,remove
#include <dirent.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "unix_cp.h"

#define BUFSIZE 4096
#define PATHSIZE 512

//--------------------------------------------------------------------
static void err_prt(const char *s1, const char *s2)
//--------------------------------------------------------------------
{
   fprintf(stderr,"unix_cp error: %s ", s1);
   perror(s2);
   exit(EXIT_FAILURE);
}

//--------------------------------------------------------------------
void copy_file(const char *file_in, const char* file_out)
//--------------------------------------------------------------------
{
   char buf[BUFSIZE];

   // open file_in
   FILE* fp_in = fopen(file_in,"r");
   if( !fp_in ) {
      err_prt("can not open ", file_in);
   }

   // create new file_out
   FILE* fp_out = fopen(file_out,"w");
   if( !fp_out ) {
      err_prt("can not create ", file_out);
   }

   // copy
   size_t nchars;
   while( (nchars = fread(buf,1, BUFSIZE,fp_in)) > 0 ) {
      if( nchars != fwrite(buf,1,nchars,fp_out) ) {
         err_prt("write to ", file_out);
      }
      if( nchars != BUFSIZE ) {
         if( feof(fp_in) ) {
            break;
         }
         err_prt("read from ", file_in);
      }
   }

   // close files
   if( fclose(fp_in) ) {
      err_prt("can not close ", file_in);
   }
   if( fclose(fp_out) ) {
      err_prt("can not close ", file_out);
   }
}

//--------------------------------------------------------------------
void copy_dir(const char *dir_in, const char* dir_out)
//--------------------------------------------------------------------
{
   char old_file[PATHSIZE]; // path to file
   char new_file[PATHSIZE];

   // open dir_in
   DIR* dp_in = opendir(dir_in);
   if( !dp_in ) {
      err_prt("can not open dir ", dir_in);
   }

   // check that dir_out exists
   struct stat statbuf;
   if( stat(dir_out,&statbuf) == 0 && S_ISDIR(statbuf.st_mode) ) {
      // printf("%s exists\n",dir_out);
    } else {
       err_prt("directory does not exist ",dir_out);
    }

   // copy only regular files
   struct dirent* de;
   while( (de = readdir(dp_in)) != NULL ) {
      // get file information
      snprintf(old_file,sizeof(old_file),"%s/%s",dir_in,de->d_name);
      if( stat(old_file,&statbuf) == -1 ) {
         continue;
      }

      if( S_ISREG(statbuf.st_mode) ) { // man stat.h
         snprintf(new_file,sizeof(new_file),"%s/%s",
               dir_out, de->d_name);
         // printf(" %s -> %s\n",old_file, new_file);
         copy_file(old_file, new_file);
      }
   }

   // close dir
   if( closedir(dp_in) < 0 ) {
      err_prt("can not close ",dir_in);
   }
}

//--------------------------------------------------------------------
void rm_whole_dir(const char *dir)
//--------------------------------------------------------------------
{
   char file[PATHSIZE];

   // open dir
   DIR* dp = opendir(dir);
   if( !dp ) {
      err_prt("can not open dir ", dir);
   }

   // remove files
   struct dirent* de;
   while( (de = readdir(dp)) != NULL ) {
      snprintf(file,sizeof(file),"%s/%s",dir,de->d_name);
      struct stat statbuf;
      if( stat(file,&statbuf) == -1 ) continue;

      // assume that all files are regular
      if( S_ISREG(statbuf.st_mode) ) { // man stat.h
         if( remove(file) ) {
            err_prt("can not remove ",file);
         }
      }
   }

   // close dir and remove it
   if( closedir(dp) < 0 ) {
      err_prt("can not close ",dir);
   }
   if( remove(dir) ) {
      err_prt("can not remove ",dir);
   }
}

//--------------------------------------------------------------------
std::string copy_dir_temp(std::string dir_in, std::string basedir)
//--------------------------------------------------------------------
{
   std::string tmp_dir = basedir + "/bean_XXXXXX";
   if ( mkdtemp((char*)tmp_dir.data()) == NULL ) {
      err_prt("can not create dir ",tmp_dir.data());
   }
   copy_dir(dir_in.data(), tmp_dir.data());

   return tmp_dir;
}

#ifdef TEST_UNIX_CP
int main(int argc, const char* argv[])
{
   // std::string dir_bd("/home/nefedov/study/bes3/bean_67/"
   std::string dir_bd("/Users/nefedov/MyHome/study/bes3/bean_67/"
         "Analysis/DatabaseSvc/dat");
   std::string basedir("/tmp");

   if ( argc==2 && atoi(argv[1])==1 ) {
      std::string dir_new = copy_dir_temp( dir_bd,basedir );
      printf("copy done:\n%s\n%s\n",dir_bd.data(),dir_new.data());
   }

   if ( argc==2 && atoi(argv[1])==2 ) {
      const char* dir_tmp = "/tmp/bean_M6pHIn";
      rm_whole_dir( dir_tmp );
      printf("remove directory %s\n",dir_tmp);
   }
}
#endif
