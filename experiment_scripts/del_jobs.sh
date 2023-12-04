vals=$(seq 2 1 9)
for v in $vals
do
echo 5980600${v}
qdel 5980600${v}
done
