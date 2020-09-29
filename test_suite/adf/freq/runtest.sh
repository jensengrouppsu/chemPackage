#! /bin/sh

echo 'Coordinates' > freq.new
echo '' >> freq.new
datadisplay.py -c freq.out >> freq.new
echo '' >> freq.new
echo '' >> freq.new
echo 'Frequencies' >> freq.new
echo '' >> freq.new
datadisplay.py freq.out >> freq.new
echo '' >> freq.new

diff freq.new freq.old > total_diffs.txt
