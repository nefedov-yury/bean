#!/bin/bash
#
# Using:
#    add_new_RootEventData.sh  path_to_boss  BOSS_version
#
# The script will create the RootEventData_BOSS_version directory and
# copy the necessary files into it.
#

# find BOSS version of RootEventData
boss_RootEventData=$(find "$1/Event/RootEventData/"\
  -maxdepth 1 -name "RootEventData*" 2>/dev/null | tail -1)
if [ ! -d "${boss_RootEventData}" ]; then
  printf "check path to BOSS: '%s'\nboss_RootEventData='%s'\n" \
    "$1" "${boss_RootEventData}"
  exit
else
  printf "RootEventData directory found in Boss: '%s'\n"\
    "${boss_RootEventData##*/}"
fi

# check destination path
bean_RootEventData="./RootEventData_$2"
if [ -d "${bean_RootEventData}" ]; then
  printf "The '%s' directory already exist, stop\n"\
    "${bean_RootEventData}"
  exit
else
  printf "The destination directory is '%s'\n"\
    "${bean_RootEventData}"
fi

# copy usin rsync (add --dry-run for test)
mkdir "${bean_RootEventData}"
rsync -av\
  --include="ChangeLog"\
  --exclude="*_rootcint.*"\
  --include="RootEventData/" --include="*.h"\
  --include="src/" --include="*.cxx"\
  --exclude="*"\
  "${boss_RootEventData}/" "${bean_RootEventData}"

# copy PROOF-INF:
# cp -r ../RootEventData_6.5.5/PROOF-INF .

# add CMakeLists.txt
cat > "${bean_RootEventData}/CMakeLists.txt" << EOF
INCLUDE( RootEventData )
EOF

printf " ... done\n"

