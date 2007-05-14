#!/bin/sh

#
# This script is used to make a new release of GETM.
# The release can be both a 'stable' and 'devel' release.
# The script should not be executed directly - but via the 
# make file targets <devel|stable>

#
# $Id
#

[ "$USER" = "kbk" ] || { echo "Only kbk can make new releases" ; exit 1; }
[ `hostname` = gate ] || { echo "Releases should be done on gate" ; exit 1; }

release_type=$1
release_version=$2

release_name=getm-$release_version
base_dir=/public/ftp/pub/getm-releases
release_dir=$base_dir/$release_type
tarfile=$release_name.tar.gz

TAG=v`echo $release_version | tr . _`
BRANCH=$TAG

RHOST=gate
RUSER=kbk
RDIR=bolding-burchard.com/src

export CVSROOT=$USER@gate:/public/cvs
export CVS_RSH=ssh

if [ -d $release_dir/$release_name ] ; then

   echo
   echo $release_name" has already been released"
   echo "update VERSION in Makefile"
   echo
   exit 1

fi

if [ "$release_type" = "stable" ] ; then
   cvs tag -b $TAG
   CVS2CL="cvs2cl -b -F $BRANCH --no-ancestors"
fi

if [ "$release_type" = "devel" ] ; then
   cvs tag $TAG
   CVS2CL="cvs2cl -F trunk"
fi

if [ "$release_type" = "branch" ] ; then
   cvs tag -b $TAG
   echo "now check out the new branch and update the Makefile"
   echo "the CVS2CL has to be modified"
   exit 0
fi


$CVS2CL && mkdir -p $release_dir/$release_name/include/ &&  mv ChangeLog VERSION $release_dir/$release_name && mv include/version.h $release_dir/$release_name/include/

cd $release_dir && cvs export -r $TAG -d $release_name getm-src && tar -cvzf $tarfile $release_name/ && rm getm-$release_type.tar.gz getm-$release_type && ln -sf $release_name.tar.gz getm-$release_type.tar.gz && ln -s $release_name getm-$release_type

scp -p $release_dir/$tarfile $RHOST:$RDIR/$release_type/
ssh $RHOST \( cd $RDIR/$release_type \; rm getm-$release_type.tar.gz \; ln -s $tarfile getm-$release_type.tar.gz \; tar -xvzf $tarfile \; rm getm-$release_type \; ln -s $release_name getm-$release_type \)

exit 0
