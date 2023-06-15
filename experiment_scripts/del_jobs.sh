vals=$(seq 1 1 6)
for v in $vals
do
echo 5669602${v}
qdel 5669602${v}
done
