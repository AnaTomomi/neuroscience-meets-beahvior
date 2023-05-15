 for n in $(cat scopus_list_as_xls2.txt);do 
	 str=$(grep $n ../tmp/scopus/*|sed 's/..\/results\/scopus\///g'|cut -d\: -f1|sed 's/%22_%22/ /g'|sed 's/%20/-/g'|sed 's/%22_[0-9]/ /g'|sed 's/[0-9]//g'|sed 's/scopus%//g'|tr ' ' '\n'|sort|uniq|tr '\n' ' ');
	 echo $n","$str 
 done 	 > ../results/scopus_keywords_as_xls.csv;

