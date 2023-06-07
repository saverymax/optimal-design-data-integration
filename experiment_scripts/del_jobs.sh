vals=$(seq 1 1 6)
for v in $vals
do
echo 5652654${v}
qdel 5652654${v}
done
