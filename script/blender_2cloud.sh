#! /bin/bash
#


readonly PROGNAME=$(basename $0)
readonly PROGDIR=$(readlink -m $(dirname $0))
readonly ARGS="$@"


usages () {
	cat <<- EOF
	usage: $PROGNAME output_directory file1 file2 ...
	Convert all input file to .qc file in output directory

	Examples:
		Convert all wrl in a directory:
		$PROGNAME output/path/ input/path/*.wrl
	EOF
}


good_args_number() {
	local nrArgs=$1
	[[ $nrArgs -ge 2 ]]
}


is_dir() {
	local path=$1
	[[ -d $path ]]
}


is_file() {
	local path=$1
	[[ -f $path ]]
}


exit_on_error() {
	local errorMsg=$1
	local errorCode=$2
	usages
	echo $errorMsg
	exit $errorCode
}


convert() {
	local in_file=$1
	local out_dir=$2
	local filename=$(basename "$in_file")
	local filename_noext="${filename%.*}"

	blender $PROGDIR/blender/empty.blend --background -noaudio --python $PROGDIR/blender_2cloud.py\
	       	-- $in_file "$out_dir/$filename_noext.qc"
}


main() {
	good_args_number $# \
		|| exit_on_error "Bad arguments numbers" 1

	local out_path=$1; shift
	local files=$@

	is_dir $out_path/ \
		|| exit_on_error "$out_path is not a directory" 2

	for file in $files
	do
		is_file $file \
			|| exit_on_error "$file is not a file" 3
	done

	for file in $files
	do
		convert $file $out_path
	done
}
main $ARGS
