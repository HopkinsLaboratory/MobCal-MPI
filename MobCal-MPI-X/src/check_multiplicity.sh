#!/bin/bash

printf "%s" $1 | grep -E "(cat|an)ion"
