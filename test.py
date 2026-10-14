
from SunStatistics.SolarCalc import System# as SC
from time import time

TX = [
    ['2026-06-21 00:00:00+00:00','2026-06-22 00:00:00+00:00','15min'],
    ['2026-06-21 00:00:00','2026-06-22 00:00:00','15min'],
    ['2026-06-21 00:00:00+00:00','2027-06-22 00:00:00+00:00','15min'],
    ['2026-06-21 00:00:00','2027-06-22 00:00:00','15min']
]
for tx in TX:
    T1 = time()
    sa=System(
        # Timestamp=['2026-06-21 00:00:00+00:00','2026-06-22 00:00:00+00:00','15min'],
        Timestamp=tx,#['2026-06-21 00:00:00','2026-06-22 00:00:00','15min'],
        Latitude=69,
        Longitude=-135,
        transmissivity=0.9,
        UTC=False,
        Slope=90,
        Aspect=0,
        # Slope = 90,
        Full_Load_Watts=24,
        Discharge_Limit=0.55,
        Batt_Sequences_in_Parallel=12
        )
    T2 = time()
    print(T2-T1)