Volconv - geometry-aware DICOM-to-Nifti volume converter
========================================================

Author: Mark White <mark@celos.net>

Source: <https://bitbucket.org/mjwhite/volconv>

Caveat
------

Always check the geometry and laterality of your images before and
after conversion; and beware that the DICOM reader is only partial.
It doesn't handle compressed transfer syntaxes at all, makes some
simplifications about reading sequences, and will probably fall over
horribly on DICOM files that aren't basically MR or CT volumes.

In particular, note that volconv has no validation for medical use.

Please feel free to send bug reports or feature requests, along with
an anonymized example of input DICOM data if that's relevant.

Installation
------------

Installing volconv shouldn't be hard.  At a minimum, you will need to
have two packages installed:

  * Python 2.5 or higher (but 3.x isn't yet supported)
  * Numpy 1.1.1 or higher

In principle, Python 2.4 and Numpy 1.0.1 are still supported, but
volconv isn't developed or tested against these any more.  It's
routinely run against Python 2.5 and 2.7, and Numpy 1.1.1 and 1.5.1.

Volconv uses internal DICOM and NIfTI functions; you don't need to
install any extra Python libraries, just Python and Numpy.

It's designed to run from the source tree.  Just unpack the archive:

    tar -xzf volconv-<version>.tar.gz

Then check it:

    cd volconv-<version>
    ./volconv

If that works (ie produces some help), you can just use it from there.

For an installation, move the directory somewhere more permanent (eg
/opt/volconv) and either add it to your path or (as root) symlink the
"volconv" executable somewhere like /usr/local/bin:

    cd /usr/local/bin
    ln -s /opt/volconv/volconv volconv

There's no reason I know of that volconv shouldn't work on Windows,
but it's not supported or tested at all.

Basics
------

Running without arguments (or with "-h" or "--help") will print usage.
This is always the most up-to-date source of documentation.  The
document you're reading now doesn't cover every feature, and is only
updated intermittently; the help is always complete, if a bit terse.

volconv is designed as a point-and-shoot tool -- it's intended to
figure out what a DICOM tree contains and write sensible Nifti files.
Typical simplest usage would be:

    volconv .

which will recursively check every file below the current directory
for DICOM header, assemble all the studies and series it finds, print
a concise list of them, then write them in the current directory with
names based on the series number and sequence name, along with an 
index file named index.json.
    
    volconv dcm -o nifti/ -p '^IM'

will recursively find all files under subdirectory 'dcm' starting with
'IM', create a new directory 'nifti' and write sensibly named volumes
and an index file there.  (Omitting the trailing slash from nifti
would cause it to be used as a prefix, rather than creating a
subdirectory.)  

Just looking
------------

Running

    volconv . -n

will report any DICOM reading problems, print the list of series, then
exit; it's useful for figuring out what's in a pile of numbered files.
(Orientations are the plane to which the image appears closest.)

    volconv . -j

will print the JSON-formatted index (useful for machine parsing).
(There's also a "-x" option to print an XML index, but this is
deprecated and contains less information than the JSON index.)

If you wish to call volconv from a script with some automation based
on the volumes found, please read this information from index.json
rather than parsing the summary written to stdout: using index.json is
both much easier and much more stable, and will give more information
(including echo times, diffusion parameters, and so on).  Most modern
languages have libraries for reading JSON -- see www.json.org for a
list, or www.mathworks.com/matlabcentral/fileexchange/23393 for a free
Matlab implementation.

    volconv path/to/file.dcm -D

will dump the DICOM header fields of a single file in a fairly
human-friendly format.  You might want to pipe through "less".  There
are several options to extract Siemens CSA private headers with -D:

    --csa-protocol
    --csa-series
    --csa-image

And some modifiers:

    -u/--dump-unknown (prints all private dicom fiels)
    --csa-full (don't omit longer CSA fields)
    --no-trunc (doesn't cut off lines at terminal width)

Selecting files
---------------

The -p/--pattern argument selects files based on name.  Saying:

    volconv . -p '^IM'

will look for files under the current directory starting with 'IM'
(the argument to -p is a regular expression; eg '.dcm$' would look for
files with the .dcm extension).

There are some other ways to filter.  You can use shell-style globs
directly like this:

    volconv dcm --fnmatch='IM*'

which is roughly equivalent to:

    volconv dcm/IM*
    
but won't fall over if there are too many matching files for your shell's
argument list.  Paths are matched relative to volconv's working directory and
may be multi-component (so "--fnmatch='series*/*.dcm'" will work); you can also
match relative to each specified starting point with --fnmatch-relative.

Often more usefully, you can use the series description (which appears
in square brackes [] in the summary using the -i/--include and
-e/--exclude options:

    volconv . -i EPI -e 2D

Will ignore any files with descriptions that don't contain EPI and any
that do contain 2D.  Exclusion always wins (so if a description
matches both -i and -e, it's excluded); given sequences:

    EPI_2D
    EPI_3D
    TSE_2D
    SCOUT

the command above would only convert "EPI_3D".  You can only specify
-e and -i once each, but the arguments are regular expressions, so
(for example) something like "-e 'SCOUT|TSE'" would exclude both.
Don't forget the single quotes when doing this from the shell!

Volconv is able to identify and handle Siemens mosaic files (sometimes
use for DTI and EPI efficiency) by looking for a label in a private
field.  If this fails, issuing -M/--csa will force the entire CSA
header to be parsed for every file, giving a definitive detection but
also slowing down the DICOM parse stage considerably.

Warnings are normally quite concise; use -v/--verbose to request more
detail, including (where relevant) the path to an example of an input
file which has generated it.  It doesn't print the name of every
individual file causing an error, because in a typical run there may
be thousands; rather, for each distinct error, the path of one input
file which has triggered it is reported.

Occasionally, you may want to convert files in the old ACR-NEMA 2.0
format (which is like an implicit DICOM without a preamble).  Volconv
will normally just reject those files; use --acr to allow it to detect
and parse them (but beware that identifying these files is a bit
flakier than DICOM, so be prepared that it may end up trying to parse
other entirely-non-DICOM files in your tree too).

Output format
-------------

By default, volconv creates Nifti files.  You can use the --gipl
option to write GIPL files (a simple internal format used by UCL's
CMIC research laboratory and its predecessors) instead.

Adding --gzip will cause output volumes to be comressed (.nii.gz).

Output files
------------

You can specify an output prefix, which might be just something to go
on the front of the filenames:

    volconv . -o phantom-

(creating files like phantom-0038-0010-0001.nii), or a directory:

    volconv . -o phantom/

which will put all the output files in a directory named "phantom",
creating it first if it doesn't exist.  Note the trailing slash:
that's important, or it's treated as a prefix.  (You can combine them
into prefixes like "phantom/test1-" to get both effects.)

When volconv encounters a volume of mixed orientation, eg a localizer
or an MPR with an orthogonal reference slice, the default behaviour is
to separate these into subseries with suffixes like sag, cor, or axi.
Including the -X/--mixed option will prevent this, concatenating the
slices but potentially creating volumes with nonsensical geometry

There are several options for output filenames.  The default, which
you can also get by specifying:

    volconv . -c

gives files with names like [seriesnumber]-[description].nii,
along with the study number (on the front) and the timepoint and echo
number (on the end) if there is more than one to distinguish.

    volconv . -t

which writes writes [seriesnumber]-[description]-[type].nii, taking
image-specific modality type from value 3 of field (0008,0008); this
might be 'M' and 'P' for magnitude and phase, for example.  You can
always find the full type field in index.json.  Finally:

    volconv . -d

is the same as -t but writes the full date and study identifier on the
front of every file, rather than adding a disambiguating series count
on the front as needed.  This makes for longer filenames, but is
useful if you want to convert several studies together; it's the only
output type which should produce consistent filenames for each
study regardless of how many studies are processed together.

There are a other less useful modes (-s, -m) described in the help.
The most flexible option is -f, which allows you to specify the naming
template yourself:

  volconv . -f 'myfile-%(ser)-%(desc)?(-t)'

will write files like myfile-0005-t2_weighted.nii, or optionally
myfile-0007-epi-0009.nii if there are multiple time points.  There's a
full list of valid substitutions in the help.  Note that you shouldn't
give the extension (.nii) as part of the the template.

If you wish to anonymize your index file, say:

    volconv . -a mystudy

If you issue -a, example DICOM filenames (ie the name of one of the
DICOM files contributing to each volume) are not written to the
index.json file.  To force them, use -E/--force-exdcm; and to force
the full path, use -F/--full-exdcm (in many offline processing
contexts, the filename or the full path will co-incidentally reveal
the patient name).

No identifiable information is written into the actual Nifti file in
any event (other than, perhaps, as part of the filename).

Match files
-----------

There's another system for simultaneously selecting which sequences to
process and defining their output names, using a "match file" to
define the relationships.  This is useful for multi-subject studies
where you wish to give short, uniform names to output files, even when
the input files may vary (for example, series numbers will be offset
if an extra localizer is run).  You use it like this:

  volconv . -w match_file.ini

The fine example_match.ini shows all the available options for this.
This is a fairly new feature which definitely still has a few bugs.

See also --symlink for writing every series then using the match
definitions to create uniform symlinks to certain volumes.

Geometry
--------

volconv attempts to do "The Right Thing" here, by default writing
Nifti images with Q-form representations corresponding to the DICOM
orientation and offset fields.

It's also possible to write Analyze-7.5-style ("A-form", with both the
Q-form and S-form codes set to zero); I don't recommend this unless
you need them, because geometry is far more ambiguous.  The option to
switch orientation representation is -O; say -Oa for A-form, -Oq for
Q-form (the default anyway).

You can flip the horizontally or vertically with -H and -V
respectively; or you can use "-r sw" or similar to move the origin to
a specific corner, assuming it starts nw in DICOM.  Flips are
accumulated, so "-r se -H" is equivalent of a single vertical flip
(h-flips cancel out).  AXIS LABELS ARE PRESERVED ON THE DATA, so with
Q-form output and a correct reader, flips will never swap patient
sides; they may not affect the display at all.  They're just altering
which pixel is transmitted first in the data array.  A-form doesn't
store orientation, so labels won't be preserved there; and remember 3D
handedness errors are possible even with correct labels.

If volconv has found a closest radiology orientation, you can say:

    volconv . -A

to reslice all the volumes approximately into the axial plane (the
closest equivalent by 90-degree rotation; it just swaps axes, no
resampling).  This may be useful for using sagittal or coronal images
with software that doesn't pay enough attention to the Q-form
orientation data.

Through-plane voxel size is taken from the gap between the first two
slices.  You can opt instead to use the recorded slice thickness
by saying:

    volconv . -g

By default, volconv uses the "Slice Location" field to stack slices in
the right order.  I've seen a few example images where this doesn't
work: in those cases, --slice-3d may solve the problem by deducing
each slice location from the "Image Position (Patient)" field.

Try also --slice-inst and --stack-unk for a couple of nastier hacks to
force slices to stack when the DICOM geometry is either broken or too
far from volconv's assumptions of regular gridding to work.

Index files
-----------

The index.json mentioned above records some parameters which would
otherwise be lost in the conversion from DICOM, including TR, TE, flip
angle, matrix size, date and time, image type and description, and
(for Siemens MR images) DTI gradient directions and SAR levels.  It's
fairly self-explanatory.

If you want SAR written in your index.json, you need to ask at
conversion time (--sar); it takes a little extra time because the full
Siemens CSA header must be parsed for each series.

One index.json file is written per conversion run containing data for
all processed sequences; it'll be over-written if you run a different
conversion with the same output prefix (-o).
