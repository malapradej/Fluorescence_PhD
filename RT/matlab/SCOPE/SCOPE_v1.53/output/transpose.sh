# A script that will transpose a tab separated file and also remove header lines.
# The script also replaces the mac (/r) line endings with unix (/n).
# <arguments> $1: the filename of the tab separated file. $2 the output filename.
# $3 is the number of header lines to remove.

mac2unix $1
if [ "$#" -eq 3 ]
then sed '1,'$3'd' $1 > temp.txt 
else cat $1 > temp.txt 
fi
awk '
{ 
    for (i=1; i<=NF; i++)  {
        a[NR,i] = $i
    }
}
NF>p { p = NF }
END {    
    for(j=1; j<=p; j++) {
        str=a[1,j]
        for(i=2; i<=NR; i++){
            str=str" "a[i,j];
        }
        print str
    }
}' temp.txt > $2
rm temp.txt

