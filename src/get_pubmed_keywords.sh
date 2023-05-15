 for n in $(cat pubmed_list_as_xls2.txt);do 
	 str=$(grep $n ../tmp/pubmeddata/*|sed 's/..\/results\/pubmeddata\///g'|cut -d\: -f1|sed 's/-/ /g'|sed 's/.txt//g'|tr ' ' '\n'|sort|uniq|tr '\n' ' ');
	 echo $n","$str 
 done 	 > ../results/pubmed_keywords_as_xls2.csv;

