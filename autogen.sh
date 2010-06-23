#!/bin/sh

# This creates all the autoconf stuff that you won't have in
# a fresh checkout of the code from the git repository

autoheader
aclocal -I m4
autoconf
automake --add-missing
