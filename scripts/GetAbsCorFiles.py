#!/usr/bin/python3
#
# This program is an analogue of the EmcShEnCalibSvc service from the
# BOSS program.
#
# PrintCalibFiles(run): Output paths to calibration files for the run,
#                       as well as the corresponding version of BOSS.
#
# GetSftParVer(run): Display SftVer, ParVer and range of 'eligible'
#                    runs for a given run and bossver (BossRelease)
# Note: range of runs are `strange'

import MySQLdb

class EmcCalibInfo:
    # ----------------------------------------------------------------
    def __init__(self):
    # ----------------------------------------------------------------
        self.verbose = False    # default: silent
        self.bossver = "7.0.9"
        self.host = "bes3db2.ihep.ac.cn"
        self.offlinedb = MySQLdb.connect (host =  self.host,
                                          user = "guest",
                                          passwd = "guestpass",
                                          db = "offlinedb")

    # Print the paths to the calibration files 'Pi0CalibFile' and
    # 'SingleGammaCalibFile' that can be used in AbsCor by calling the
    # ReadDatac3p('Pi0CalibFile') and
    # ReadCorFunpara('SingleGammaCalibFile') functions, respectively.
    # ----------------------------------------------------------------
    def PrintCalibFiles(self,run):
    # ----------------------------------------------------------------
        cur_offline = self.offlinedb.cursor()
        sql = "SELECT SftVer,singleGammaCalib,pi0Calib"\
              " FROM EmcShEnCalibConst"\
              " WHERE RunFrom <= '{RunFrom}' and RunTo >= '{RunTo}'"
        sql = sql.format(RunFrom=run,RunTo=run)
        if self.verbose:
            print(sql)
        cur_offline.execute (sql)
        if self.verbose:
            print(" count= ",cur_offline.rowcount)
        print("run#",run,":")
        if cur_offline.rowcount == 0:
            print("ERROR: can not find records")
        else:
            so = "     BossVer= {:15s}"\
                 " (Pi0CalibFile,SingleGammaCalibFile)\n"\
                 "     {}\n"\
                 "     {}\n"
            for row in cur_offline.fetchall():
                (SftVer,SingleGammaCalibFile,Pi0CalibFile)=row
                print(so.format(SftVer,Pi0CalibFile,SingleGammaCalibFile))
        cur_offline.close()

    # Get SftVer, ParVer and range of 'eligible' runs
    # for a given run and bossver (BossRelease)
    # ----------------------------------------------------------------
    def GetSftParVer(self,run):
    # ----------------------------------------------------------------
        cur_offline = self.offlinedb.cursor()
        sql = "SELECT RunFrom,RunTo,SftVer,ParVer"\
              " FROM CalVtxLumVer"\
              " WHERE BossRelease='{BossRelease}' and"\
              " RunFrom <= '{RunFrom}' and RunTo >= '{RunTo}'"\
              " and DataType='EmcShEnCalib'"
        sql = sql.format(BossRelease=self.bossver,
                         RunFrom=run,RunTo=run)
        if self.verbose:
            print(sql)
        cur_offline.execute (sql)
        if cur_offline.rowcount == 1:
            (RunFrom,RunTo,SftVer,ParVer) = cur_offline.fetchone()
            os="runs: [{},{}] -> BossVer= {} ParVer= {}"
            print(os.format(RunFrom,RunTo,SftVer,ParVer))
        else:
            print("ERROR: can not find records\n",sql)
        cur_offline.close()

# Examples to use ----------------------------------------------------
r = EmcCalibInfo()
#  r.verbose = True
r.bossver = "7.0.9"

# Psi(2S) 2009 : run8093-run9025
#  r.PrintCalibFiles(8093)

# Psi(2S) 2012: run25338-run27090
#  r.PrintCalibFiles(25338)

# Psi(2S) 2021: run66257-run69292
r.PrintCalibFiles(66257)
#  r.GetSftParVer(66257)

# 3650 2009: run 9613 - 9779
# 3650 2012: run 33725 - 33772
# 3650 2022: run 69612 - 70132
#  r.PrintCalibFiles(9613)
#  r.PrintCalibFiles(33725)
# r.PrintCalibFiles(69612)
