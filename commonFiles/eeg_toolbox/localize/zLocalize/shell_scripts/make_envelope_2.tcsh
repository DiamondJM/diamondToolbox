#! /bin/tcsh -f
#
# make_envelope.tcsh
#
# An adaption by John C. Cocjin on 29 July 2016 which is a subset of the script
# mris_comput_lgi, created originally by Marie Schaer
#
#
# Computes local measurements of gyrification at points over cortical surface.
# --help option will show usage
#
# This script implements the work of Marie Schaer et al., as described in:
#
# "A Surface-based Approach to Quantify Local Cortical Gyrification",
# Schaer M. et al., IEEE Transactions on Medical Imaging, 2007, TMI-2007-0180
#
# Original Author: Marie Schaer
# Original filename: mris_compute_lgi
#
# Copyright (C) 2007-2008,
# The General Hospital Corporation (Boston, MA).
# All rights reserved.
#
# Distribution, usage and copying of this software is covered under the
# terms found in the License Agreement file named 'COPYING' found in the
# FreeSurfer source code root directory, and duplicated here:
#FreeSurfer Software License Agreement ("Agreement")
#Version 1.0 (February 2011)
#
#This Agreement covers contributions to and downloads from the
#FreeSurfer project ("FreeSurfer") maintained by The General Hospital
#Corporation, Boston MA, USA ("MGH"). Part A of this Agreement applies to
#contributions of software and/or data to FreeSurfer (including making
#revisions of or additions to code and/or data already in FreeSurfer). Part
#B of this Agreement applies to downloads of software and/or data from
#FreeSurfer. Part C of this Agreement applies to all transactions with
#FreeSurfer. If you distribute Software (as defined below) downloaded from
#FreeSurfer, all of the paragraphs of Part B of this Agreement must be
#included with and apply to such Software.
#
#Your contribution of software and/or data to FreeSurfer (including prior
#to the date of the first publication of this Agreement, each a
#"Contribution") and/or downloading, copying, modifying, displaying,
#distributing or use of any software and/or data from FreeSurfer
#(collectively, the "Software") constitutes acceptance of all of the
#terms and conditions of this Agreement. If you do not agree to such
#terms and conditions, you have no right to contribute your
#Contribution, or to download, copy, modify, display, distribute or use
#the Software.
#
#PART A. CONTRIBUTION AGREEMENT - License to MGH with Right to Sublicense
#("Contribution Agreement").
#
#1. As used in this Contribution Agreement, "you" means the individual
#contributing the Contribution to FreeSurfer and the institution or
#entity which employs or is otherwise affiliated with such
#individual in connection with such Contribution.
#
#2. This Contribution Agreement applies to all Contributions made to
#FreeSurfer, including without limitation Contributions made prior to
#the date of first publication of this Agreement. If at any time you
#make a Contribution to FreeSurfer, you represent that (i) you are
#legally authorized and entitled to make such Contribution and to
#grant all licenses granted in this Contribution Agreement with
#respect to such Contribution; (ii) if your Contribution includes
#any patient data, all such data is de-identified in accordance with
#U.S. confidentiality and security laws and requirements, including
#but not limited to the Health Insurance Portability and
#Accountability Act (HIPAA) and its regulations, and your disclosure
#of such data for the purposes contemplated by this Agreement is
#properly authorized and in compliance with all applicable laws and
#regulations; and (iii) you have preserved in the Contribution all
#applicable attributions, copyright notices and licenses for any
#third party software or data included in the Contribution.
#
#3. Except for the licenses granted in this Agreement, you reserve all
#right, title and interest in your Contribution.
#
#4. You hereby grant to MGH, with the right to sublicense, a
#perpetual, worldwide, non-exclusive, no charge, royalty-free,
#irrevocable license to use, reproduce, make derivative works of,
#display and distribute the Contribution. If your Contribution is
#protected by patent, you hereby grant to MGH, with the right to
#sublicense, a perpetual, worldwide, non-exclusive, no-charge,
#royalty-free, irrevocable license under your interest in patent
#rights covering the Contribution, to make, have made, use, sell and
#otherwise transfer your Contribution, alone or in combination with
#any other code.
#
#5. You acknowledge and agree that MGH may incorporate your
#Contribution into FreeSurfer and may make FreeSurfer available to members
#of the public on an open source basis under terms substantially in
#accordance with the Software License set forth in Part B of this
#Agreement. You further acknowledge and agree that MGH shall
#have no liability arising in connection with claims resulting from
#your breach of any of the terms of this Agreement.
#
#6. YOU WARRANT THAT TO THE BEST OF YOUR KNOWLEDGE YOUR CONTRIBUTION
#DOES NOT CONTAIN ANY CODE THAT REQURES OR PRESCRIBES AN "OPEN
#SOURCE LICENSE" FOR DERIVATIVE WORKS (by way of non-limiting
#example, the GNU General Public License or other so-called
#"reciprocal" license that requires any derived work to be licensed
#under the GNU General Public License or other "open source
#license").
#
#PART B. DOWNLOADING AGREEMENT - License from MGH with Right to Sublicense
#("Software License").
#
#1. As used in this Software License, "you" means the individual
#downloading and/or using, reproducing, modifying, displaying and/or
#distributing the Software and the institution or entity which
#employs or is otherwise affiliated with such individual in
#connection therewith. The General Hospital Corporation ("MGH")
#hereby grants you, with right to sublicense, with
#respect to MGH's rights in the software, and data, if any,
#which is the subject of this Software License (collectively, the
#"Software"), a royalty-free, non-exclusive license to use,
#reproduce, make derivative works of, display and distribute the
#Software, provided that:
#
#(a) you accept and adhere to all of the terms and conditions of this
#Software License;
#
#(b) in connection with any copy of or sublicense of all or any portion
#of the Software, all of the terms and conditions in this Software
#License shall appear in and shall apply to such copy and such
#sublicense, including without limitation all source and executable
#forms and on any user documentation, prefaced with the following
#words: "All or portions of this licensed product (such portions are
#the "Software") have been obtained under license from The General Hospital
#Corporation "MGH" and are subject to the following terms and
#conditions:"
#
#(c) you preserve and maintain all applicable attributions, copyright
#notices and licenses included in or applicable to the Software;
#
#(d) modified versions of the Software must be clearly identified and
#marked as such, and must not be misrepresented as being the original
#Software; and
#
#(e) you consider making, but are under no obligation to make, the
#source code of any of your modifications to the Software freely
#available to others on an open source basis.
#
#2. The license granted in this Software License includes without
#limitation the right to (i) incorporate the Software into
#proprietary programs (subject to any restrictions applicable to
#such programs), (ii) add your own copyright statement to your
#modifications of the Software, and (iii) provide additional or
#different license terms and conditions in your sublicenses of
#modifications of the Software; provided that in each case your use,
#reproduction or distribution of such modifications otherwise
#complies with the conditions stated in this Software License.
#
#3. This Software License does not grant any rights with respect to
#third party software, except those rights that MGH has been
#authorized by a third party to grant to you, and accordingly you
#are solely responsible for (i) obtaining any permissions from third
#parties that you need to use, reproduce, make derivative works of,
#display and distribute the Software, and (ii) informing your
#sublicensees, including without limitation your end-users, of their
#obligations to secure any such required permissions.
#
#4. The Software has been designed for research purposes only and has
#not been reviewed or approved by the Food and Drug Administration
#or by any other agency. YOU ACKNOWLEDGE AND AGREE THAT CLINICAL
#APPLICATIONS ARE NEITHER RECOMMENDED NOR ADVISED. Any
#commercialization of the Software is at the sole risk of the party
#or parties engaged in such commercialization. You further agree to
#use, reproduce, make derivative works of, display and distribute
#the Software in compliance with all applicable governmental laws,
#regulations and orders, including without limitation those relating
#to export and import control.
#
#5. The Software is provided "AS IS" and neither MGH nor any
#contributor to the software (each a "Contributor") shall have any
#obligation to provide maintenance, support, updates, enhancements
#or modifications thereto. MGH AND ALL CONTRIBUTORS SPECIFICALLY
#DISCLAIM ALL EXPRESS AND IMPLIED WARRANTIES OF ANY KIND INCLUDING,
#BUT NOT LIMITED TO, ANY WARRANTIES OF MERCHANTABILITY, FITNESS FOR
#A PARTICULAR PURPOSE AND NON-INFRINGEMENT. IN NO EVENT SHALL
#MGH OR ANY CONTRIBUTOR BE LIABLE TO ANY PARTY FOR DIRECT,
#INDIRECT, SPECIAL, INCIDENTAL, EXEMPLARY OR CONSEQUENTIAL DAMAGES
#HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY ARISING IN ANY WAY
#RELATED TO THE SOFTWARE, EVEN IF MGH OR ANY CONTRIBUTOR HAS
#BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGES. TO THE MAXIMUM
#EXTENT NOT PROHIBITED BY LAW OR REGULATION, YOU FURTHER ASSUME ALL
#LIABILITY FOR YOUR USE, REPRODUCTION, MAKING OF DERIVATIVE WORKS,
#DISPLAY, LICENSE OR DISTRIBUTION OF THE SOFTWARE AND AGREE TO
#INDEMNIFY AND HOLD HARMLESS MGH AND ALL CONTRIBUTORS FROM AND
#AGAINST ANY AND ALL CLAIMS, SUITS, ACTIONS, DEMANDS AND JUDGMENTS
#ARISING THEREFROM.
#
#6. None of the names, logos or trademarks of MGH or any of
#MGH's affiliates or any of the Contributors, or any funding
#agency, may be used to endorse or promote products produced in
#whole or in part by operation of the Software or derived from or
#based on the Software without specific prior written permission
#from the applicable party.
#
#7. Any use, reproduction or distribution of the Software which is not
#in accordance with this Software License shall automatically revoke
#all rights granted to you under this Software License and render
#Paragraphs 1 and 2 of this Software License null and void.
#
#8. This Software License does not grant any rights in or to any
#intellectual property owned by MGH or any Contributor except
#those rights expressly granted hereunder.
#
#PART C. MISCELLANEOUS
#
#This Agreement shall be governed by and construed in accordance with
#the laws of The Commonwealth of Massachusetts without regard to
#principles of conflicts of law. This Agreement shall supercede and
#replace any license terms that you may have agreed to previously with
#respect to FreeSurfer.
#

set VERSION = '$Id: mris_compute_lgi,v 1.23 2009/01/26 22:14:05 nicks Exp $'
set PrintHelp = 0;
set RunIt = 1;
set cmdargs = ($argv);
set use_mris_extract = 1
set closespheresize = 15
set smoothiters = 30
set radius = 25
set stepsize = 100
#set echo=1
set start=`date`
if($#argv == 0) then
# zero args is not allowed
goto usage_exit;
endif
goto parse_args;
parse_args_return:
goto check_params;
check_params_return:
# begin...
#---------

# temporary work files go here...
set tmpdir = ($PWD/tmp-mris_compute_lgi-${input})

set outer =  ${tmpdir}/${input}-outer

if ( $RunIt && ! -e ${outer} ) then
echo "ERROR: make_outer_surface did not create output file '${outer}'!"
exit 1
endif
#
# mris_extract_main_component (optional, default)
#
if ($use_mris_extract) then
set cmd=(mris_extract_main_component \
${tmpdir}/${input}-outer \
${tmpdir}/${input}-outer-main)
echo "================="
echo "$cmd"
echo "================="
if ($RunIt) $cmd
if($status) then
echo "ERROR: $cmd failed!"
exit 1;
endif
else
set cmd=(cp ${tmpdir}/${input}-outer ${tmpdir}/${input}-outer-main)
echo "================="
echo "$cmd"
echo "================="
if ($RunIt) $cmd
if($status) then
echo "ERROR: $cmd failed!"
exit 1;
endif
endif
#
# mris_smooth
#
# smooth this jaggy, tessellated surface
# Mike adding the area flag to see if this stops erosion
echo "attemption area normalization"
set cmd=(mris_smooth -area -nw -n ${smoothiters} \
${tmpdir}/${input}-outer-main \
./${input}-outer-smoothed)
echo "================="
echo "$cmd"
echo "================="
if ($RunIt) $cmd
if($status) then
echo "ERROR: $cmd failed!"
exit 1;
endif
#
# mris_euler_number (a QA check, total defects should = 0)
#
set cmd=(mris_euler_number ./${input}-outer-smoothed)
echo "================="
echo "$cmd"
echo "================="
if ($RunIt) $cmd
if($status) then
echo "ERROR: $cmd failed!"
exit 1;
endif
#
# mris_convert
#
# output normals of the smoothed outer surface
set cmd=(mris_convert -n \
${input}-outer-smoothed \
${tmpdir}/${input}-outer-smoothed-normals.asc)
echo "================="
echo "$cmd"
echo "================="
if ($RunIt) $cmd
if($status) then
echo "ERROR: $cmd failed!"
exit 1;
endif
#--------
# end...
# *MST 02/17 remove deletion of tmp directory
#set cmd=(rm -Rf $tmpdir)
#echo "================="
#echo "$cmd"
#echo "================="
#if ($RunIt) $cmd
set end=`date`
echo "done."
echo "Start: $start"
echo "End:   $end"
exit 0
#----------------------------------------------------------#
############--------------##################
parse_args:
set cmdline = ($argv);
while( $#argv != 0 )
set flag = $argv[1]; shift;
switch($flag)
case "--close_sphere_size":
set closespheresize = $argv[1]; shift;
#echo "Using sphere size of ${closespheresize}mm for morph closing op."
breaksw
case "--smooth_iters":
set smoothiters = $argv[1]; shift;
echo "Smoothing using ${smoothiters} iterations."
breaksw
case "--step_size":
set stepsize = $argv[1]; shift;
echo "Skipping every ${stepsize} vertices during lGI calcs."
breaksw
case "--help":
set PrintHelp = 1;
goto usage_exit;
exit 0;
breaksw
case "--version":
echo $VERSION
exit 0;
breaksw
case "-i":
case "--i":
case "--input":
if ( $#argv == 0) goto arg1err;
set input = $argv[1]; shift;
#echo ${input}
breaksw
case "--dont_extract":
set use_mris_extract = 0
breaksw
case "--debug":
case "--echo":
set echo = 1;
set verbose = 1
breaksw
case "--dontrun":
set RunIt = 0;
breaksw
default:
breaksw
endsw
end
goto parse_args_return;
############--------------##################
############--------------##################
arg1err:
echo "ERROR: flag $flag requires one argument"
exit 1
############--------------##################
############--------------##################
check_params:
if(! $?FREESURFER_HOME ) then
echo "ERROR: environment variable FREESURFER_HOME not set."
exit 1;
endif
if(! -e $FREESURFER_HOME ) then
echo "ERROR: FREESURFER_HOME $FREESURFER_HOME does not exist."
exit 1;
endif
if(! $?input) then
echo "ERROR: missing input filename!  See  mris_compute_lgi --help"
exit 1;
endif
if(! -e ${input} ) then
echo "ERROR: input file '${input}' does not exist."
exit 1;
endif
goto check_params_return;
############--------------##################
############--------------##################
usage_exit:
echo ""
echo "USAGE: mris_compute_lgi [options] --i <input surface>"
echo ""
echo "Produces a surface map file containing local gyrification measures."
echo "Output file is named <input surface>_lgi, where <input surface> is the"
echo "specified input surface (ex. lh.pial produces lh.pial_lgi)."
echo ""
echo "Required Arguments"
echo "  --i       : input surface file, typically lh.pial or rh.pial"
echo ""
echo "Optional Arguments"
echo "  --close_sphere_size <mm> : use sphere of size <mm> mm for morph"
echo "                             closing operation (default: ${closespheresize})"
echo "  --smooth_iters <iters>   : smooth outer-surface <iters> number of"
echo "                             iterations (default: ${smoothiters})"
echo "  --step_size <steps>      : skip every <steps> vertices when"
echo "                             computing lGI (default: ${stepsize})"
echo "  --help    : short descriptive help"
echo "  --version : script version info"
echo "  --echo    : enable command echo, for debug"
echo "  --debug   : same as --echo"
echo "  --dontrun : just show commands (dont run them)"
echo ""
if(! $PrintHelp) exit 1;
echo Version: $VERSION
cat $0 | awk 'BEGIN{prt=0}{if(prt) print $0; if($1 == "BEGINHELP") prt = 1 }'
exit 1;
#---- Everything below here is printed out as part of help -----#
BEGINHELP
Computes local measurements of gyrification at thousands of points over the
entire cortical surface using the method described in:
"A Surface-based Approach to Quantify Local Cortical Gyrification",
Schaer M. et al., IEEE Transactions on Medical Imaging, 2007, TMI-2007-0180
Input is a pial surface mesh, and the output a scalar data file containing
the local gyrification index data at each vertices.
Example:
mris_compute_lgi --i lh.pial
produces lh.pial_lgi
See also http://surfer.nmr.mgh.harvard.edu/fswiki/LGI

