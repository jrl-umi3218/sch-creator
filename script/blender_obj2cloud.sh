#! /bin/sh

local_dir=`dirname $0`

if [ $# -lt 2 ]
then
  echo "Usage: `basename $0` input_path output_path"
  exit 1
fi

in_path=$1
out_path=$2

if [ ! -d "$in_path" ]
then
  echo "'$in_path' is not a valid directory."
  exit 2
fi

echo "$out_path"
if [ ! -d "$out_path" ]
then
  echo "'$out_path' is not a valid directory."
  exit 2
fi

for file in $in_path/*.obj
do
  echo $file
  file_name=`basename $file .obj`

  blender $local_dir/blender/empty.blend --background -noaudio --python $local_dir/blender_obj2cloud.py --  $file "$out_path/$file_name.qc"

done

