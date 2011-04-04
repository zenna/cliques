FILES=/home/zenna/repos/uniproj/data/Aaron/*
for f in $FILES
do
  #echo "Processing $f file..."
  base=$(basename $f)
  python matrix_to_matlab.py $f > /home/zenna/repos/uniproj/data/brains30/${base}
  #f="${f%.*}"
  #echo "bb $f"
  #./find_parts $f
  #counter=$(wc -l < ${f}_unsorted.stabilities )
  #sort -g ${f}_unsorted.stabilities > ${f}.stabilities
  #sed -i "1 i $counter" ${f}.stabilities
  #rm ${f}_unsorted.stabilities
  #echo “extension: ${f##*.}”
done

