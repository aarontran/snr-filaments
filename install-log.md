Ad hoc install log (June 19 2014)
=================================
(finished all installs June 19...)

In general I am installing most astro software to `~/bin/`.
This should be added to `$PATH` in `.bash_profile`.

1. Set Finder, Firefox/Chrome, Terminal preferences
2. Set System Preferences (remap Caps Lock key to Ctrl)
3. RC files -- import/reuse/twiddle the following:

    .bash_profile
    .screenrc
    .vimrc
    .vim/colors/vividchalk.vim
    .matplotlib/matplotlibrc
    .xspec/Xspec.init  # Just need to change the help doc viewers

4. Software (in approximate order of installation)

* MacVim (update `.bash_profile` alias accordingly)
* TextMate (enable shell support for `mate` in preferences)
* XQUARTZ 2.7.6
* XCode (+ command line tools).  Used aaron.tran@nasa.gov AppleID.
* CIAO + smoke tests
* homebrew
* MacTeX 2013
* `brew install git`  # Update to 2.0.0
* `brew install gcc --universal`  # 32 bit gfortran for HEASOFT

5. Brew/pip for Python.  Tried Anaconda, but there's no 32 bit support.

* `brew install python --universal --framework`
  (pip packages go into `/usr/local/lib/python2.7/site-packages`)
  (set `$PATH` in `.bash_profile` to prefer /usr/local/bin over /usr/bin)
  (python header files may be in the wrong places for PyXspec, symlink?)
  (important note -- set $PATH variable before any further package installs)

* pip install numpy
* pip install ipython[all]  # Good with `iptest` but not all packages vetted
* pip install mock  # For ipython tests

* pip install scipy  # Prerequisite: gcc / gfortran
* brew install freetype
* pip install matplotlib  # Prerequisite: freetype
* pip install pandas  # For good measure
* pip install sympy  # For good measure

Note: due to some trouble with libpng (matplotlib builds on older XQuartz
libpng? but iPython uses newer brewed version... I try to force matplotlib to
build on the new libpng, following this stackoverflow
[http://stackoverflow.com/a/23917785](answer)
(see also: [http://stackoverflow.com/q/22898245](another Q),
[https://github.com/ipython/ipython/issues/5228](github)).


Other stuff
-----------

* PyDS9 1.7 (/usr/local/lib/python2.7/site-packages)  # Prerequisite: Python
* Funtools 1.4.4 installed to `/Users/atran3/bin/saord`


HEASOFT 6.15.1 install
----------------------

Check compiler locations/versions
1. gfortran4.8.3 (from brewed gcc4.8.3)
2. perl 5.12.4 (from apple)

gcc, g++, perl are from /usr/bin
But, Python is brewed in /usr/local/bin

Add symlinks manually for Python header/library files, as follows:

    /usr/local/lib$ ln -s
    ../Cellar/python/2.7.7_2/Frameworks/Python.framework/Versions/2.7/lib/libpython2.7.dylib
    libpython2.7.dylib

    /usr/local/include$ ln -s
    ../Cellar/python/2.7.7_2/Frameworks/Python.framework/Versions/2.7/include/python2.7
    python2.7

Now execute the configure/make/make install

    ./configure --prefix=/Users/atran3/bin/heasoft-6.15.1/
    make
    make install


Updating file permissions (from FAT32 USB key...)
-------------------------------------------------

To update directories/files, to not have write/read access by default for
group/other and stuff...
    find . -type f -exec chmod 644 {} \;
    find . -type d -exec chmod 755 {} \;
Then allow user execution for shell files:
    find . -name .sh -exec chmod u+x {} \;

For reference:
6 = rw-
4 = r--
7 = rwx
5 = r-x

