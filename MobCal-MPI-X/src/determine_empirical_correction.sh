#Determine if ORCA file should have MobCal empirical corrections applied

THRESHOLD=150
if [[ -z $1 ]]; then
    printf "Input file value is a required input.\n" >&2
    exit 1
elif [[ -z $2 ]] || ! [[ $2 -gt 0 ]]; then
    printf "Threshold value not provided. Using '%s' as default. \n" $THRESHOLD >&2
else
    THRESHOLD=$2
fi

IONMASSLINE=(`cat "$1" | grep "Total Mass"`)
IONMASS=${IONMASSLINE[-2]}

printf "%s >= %s\n" $IONMASS $THRESHOLD | bc -l
