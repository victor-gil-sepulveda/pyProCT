import prody
from pyproct.driver.time.timerHandler import TimerHandler

# Set current version
__version__ = "1.7.1"

# Set prody verbosity to avoid excess of cmd line messages
prody.confProDy(verbosity="none")
