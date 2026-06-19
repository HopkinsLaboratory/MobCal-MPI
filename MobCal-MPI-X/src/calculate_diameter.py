import sys
import math

#Declare value
inputValue = 0


if (len(sys.argv) > 1):
    inputValue = float(sys.argv[1])

computation = round(((6*inputValue)/math.pi) ** (1/3), 3)    
print(computation)

