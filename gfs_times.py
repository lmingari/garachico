#!/usr/bin/env python
import argparse
import logging
import math
from datetime import datetime, timedelta 
import json
try:
    from datetime import UTC
except ImportError:
    from datetime import timezone
    UTC = timezone.utc

# Configure basic logging
logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger("GFS_TIME")

def main(args):
    # Reference date
    date = datetime.strptime(args.date,'%Y%m%d').replace(tzinfo=UTC)
    # Latency for availability of GFS data
    latency = timedelta(hours=args.latency)
    # Current datetime
    now = datetime.now(UTC)
    today = datetime(now.year, now.month, now.day, tzinfo=UTC)

    # Start datetime for searching
    t0_search = now - latency

    # Start and end datetimes for model
    t0_model = date + timedelta(hours=args.hinit)
    t1_model = date + timedelta(hours=args.hend)

    # Assume there is data available for 10 days
    foundTime = False
    data = {}
    for days in [timedelta(days=i) for i in range(10)]:
        if foundTime: break
        for cycle in [18,12,6,0]:
            t = today - days + timedelta(hours=cycle)
            if t <= t0_model and t < t0_search:
                foundTime = True
                data['cycle'] = cycle
                data['date'] = t.strftime("%Y%m%d")
                break

    if foundTime:
        logger.info(f"Found forecasts")
    else:
        raise Exception("GFS data not available for this time period")

    # Compute time range
    tmin = int(math.floor((t0_model - t).total_seconds()/3600))
    tmax = int(math.ceil ((t1_model - t).total_seconds()/3600))
    tmax += (tmax-tmin)%args.step

    # Assume forecasts for 120 hours
    if tmin < tmax and tmax <= 120:
        data['tmin'] = tmin
        data['tmax'] = tmax
    else:
        raise Exception("Request out of time range")

    print(json.dumps(data))

if __name__ == "__main__":
    #
    # Argument parser
    #
    parser=argparse.ArgumentParser(description="Compute times for GFS request")
    parser.add_argument("--date", required=True, metavar='YYYYMMDD', help='Reference date in format YYYYMMDD')
    parser.add_argument('--hinit', required=True, metavar='start_time', type=float, help='Start time in hours since date at 00:00')
    parser.add_argument('--hend', required=True, metavar='end_time', type=float, help='End time in hours since date at 00:00')
    parser.add_argument('--step', metavar='time_resolution', type=int, default = 1, help='Time step for forecast hours')
    parser.add_argument('--latency', metavar='latency', type=int, default = 6, help='Latency time for GFS availability')
    args=parser.parse_args()
    #
    # Main program
    #
    logger.info("Getting times for GFS request...")
    main(args)
    logger.info("Done!")
