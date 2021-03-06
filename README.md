# BCH_assignment

To identify hotspots, I first began by splitting chromosomes up into 2 Mb bins. I got the idea to split the regions into 2 Mb bins from the provided ChromInfo.txt file and from a paper from Frederick W. Walt in 2012: https://doi.org/10.1016/j.cell.2011.07.049. A major limitation of this approach is that the last bin of each chromosome is not exactly 2 Mb long. Also, a potential hotspot may be missed if it is split up into two insignificant bins. Nevertheless, I then found Q3 and Q1 values and obtained the IQR. I defined a hotsplot as a bin with translocatinos exceeding Q75+(1.5 x IQR). I got this idea from https://www.scribbr.com/statistics/outliers/. Also, I borrowed the ReadFile function from rebelCoderBio on youtube: https://www.youtube.com/c/rebelCoderBio, this is a funciton I use frequently, but credit goes to him nevertheless. 

The output of my script is a dataframe of hotspots containing the chromosome name, start position, end position, and number of translocations. The script also creates a graph which is useful for visualizing the number of translocations in each bin across the genome. Many bins had zero translocations, most bins had 5-20 translocations, and few bins were considered hotspots.

If you have any questions about the code, feel free to ask me at donniekramer@gmail.com
