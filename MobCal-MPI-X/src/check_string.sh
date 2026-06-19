#!/bin/bash

printf "%s" $1 | grep -E "(deprotonate|anion)"
