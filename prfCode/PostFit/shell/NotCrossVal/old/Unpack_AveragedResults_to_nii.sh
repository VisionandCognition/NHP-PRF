#/bin/bash
MONKEY=$1
AR=../Results/${MONKEY}/AveragedResults_SD
echo 'Monkey is'${MONKEY}
for TF in ${AR}/Thr_*; do
	echo 'Doing folder '%TF
	mkdir -p ${TF}/nii
	for f in $TF/*.nii.gz; do
		ff=$(basename "$f")
		fff=$(basename "$f" .gz)
		echo 'file '${ff}
		gunzip -ck $f > ${TF}/nii/${fff}
	done
done
