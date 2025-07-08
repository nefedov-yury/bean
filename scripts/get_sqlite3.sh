#!/bin/sh
# https://sqlite.org/howtocompile.html
# + sqlite3.c: The SQLite amalgamation source file
# + sqlite3.h: The header files that accompanies sqlite3.c and defines
#              the C-language interfaces to SQLite.
# - shell.c:   The command-line interface program with main() routine.
# - sqlite3ext.h: Programming Loadable Extensions
#                 (https://www.sqlite.org/loadext.html)

#  SQLITE_URL="https://www2.sqlite.org/2022/sqlite-amalgamation-3390200.zip"
SQLITE_URL="https://sqlite.org/2024/sqlite-amalgamation-3450300.zip"
LOCAL_ZIP="sqlite3-amalgamation.zip"
if [ ! -f ${LOCAL_ZIP} ]; then
  wget ${SQLITE_URL} -O ${LOCAL_ZIP}
fi
unzip -u -j ${LOCAL_ZIP} "sqlite-amalgamation-*/sqlite3.*"
