Original:

Apply lasieve5_ggnfs.patch
Apply lasieve5_boinc.patch

From freebsd build:
Apply patch-general
Apply patch-nodebug to remove debug symbols

patch-freebsd is for FreeBSD builds
patch-very_large_q might extend sieving range
patch-gcc46 makes it compile with new gcc
patch-clang makes it compile with new clang ... maybe

smjsmod - mods from Stephen Searle to enable Mac and Windows builds.
Attached is a tar file of my modified lasieve 5 source. I haven’t stripped out all the comments I’ve added to 
keep track of my changes - I’ll remove most of those at some point. I’ve modified the Makefiles a bit. 
Typing just ‘make’ should build all the ‘e’ versions (gnfs-lasieve4I11e etc). Other new targets include 
‘allf’ to build all the ‘f’ versions,  ‘allg’ to build all the ‘g’ versions. You should NOT need to build 
the libraries in the asm directory before typing make in the top directory. I’ve added a few targets for cleaning up:
   ‘realclean’ - delete all c and h intermediate build files of the web source, in addition to the standard ‘clean'
   ‘binclean’  - delete all the compiled binaries and .out files from the tests in addition to ‘realclean’.

Before you build on mac or linux you need to:
ln -s athlon64 asm
On windows it may just be easier to make an asm directory and copy the athlon64 files to it

and:
Edit the Makefile and asm/Makefile to set the compile flags you want and the compiler and ctangle commands.
Edit asm/ls-defs.asm to specify the OS type (linux, osx or windows).  l1_bits is set to 15, which seems to 
be optimal on my processors.

For mac on OSX 10.9 there is no gcc installed by default. I installed gcc 4.9 using brew and it produces 
slightly faster binaries compared to the clang compiler so you might want to try that. 

Notes:
The f and g versions haven’t been adapted for GGNFS format input, so I’ve used the original input-poly 
function for them (renamed to input_poly_orig).

I’ve added a define (NO_TD_CLOCK) to asm/siever-config.w to allow defining out all the clock calls which 
saves a bit of time. This isn’t defined by default so you have times.

The mingw compile setup was basically as described in the mersenne post by Dan Ee except I used more 
recent versions of mingw64 compiler and msys, and installed cweb to process the .w files. I built MPIR.

The windows version does not handle ctrl-c very well (it occasionally succeeds in writing the last_spq 
file but mostly fails. This happened with the version compiled by Dan Ee too. My guess is the setjmp 
isn’t able to handle jumping out of the functions which are using the sysv calling convention calls, 
but I haven’t investigated it as its a minor issue.

There are a couple of directories containing test poly files (forumexs and readmeex). 

The xeon64 directory has been adapted to work with osx as well as linux, but I haven’t tried to make 
that directory work with windows. The athlon64 directory seems like its the most recent so that’s the 
one I focused on most.

I compiled against gmp version 6 (the latest one) on osx. On windows I used MPIR. 

The CFLAGS for osx using gcc-4.9 uses ‘-Wa,-q’. That tells it to use the clang as command rather than 
the one in /usr/bin. That option probably isn’t valid on linux, but its important if you’re using osx 10.9
 
Let me know how you get on.

Regards

Steve

