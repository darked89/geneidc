#!/bin/bash

uncrustify -c ./uncrustify.cfg --replace  ../include/*.h

uncrustify -c ./uncrustify.cfg --replace  ../src/*.c
