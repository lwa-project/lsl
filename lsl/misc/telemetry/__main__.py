import argparse

from lsl.misc.telemetry import _INSTALL_KEY, enable, disable, is_active

parser = argparse.ArgumentParser(
        description='update the LSL telemetry setting',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
tgroup = parser.add_mutually_exclusive_group(required=False)
tgroup.add_argument('-e', '--enable', action='store_true', 
                    help='enable telemetry for LSL')
tgroup.add_argument('-d', '--disable', action='store_true', 
                    help='disable telemetry for LSL')
parser.add_argument('-k', '--key', action='store_true',
                    help='show install identification key')
args = parser.parse_args()

# Toggle
if args.enable:
    enable()
elif args.disable:
    disable()
    
# Report
## Status
print("LSL Telemetry is %s" % ('active' if is_active() else 'in-active'))

## Key
if args.key:
    print("  Identification key: %s" % _INSTALL_KEY)
