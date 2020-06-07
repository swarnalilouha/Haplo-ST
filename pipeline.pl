#!/bin/perl
use LWP::Simple;
use Term::ANSIColor;

print "Enter the value for the 'length criteria'\n";
chop($length_criteria=<STDIN>);
system("mkdir ./cleaned_data");
opendir(DIR1,raw_data);
@file_list1=readdir(DIR1);
closedir(DIR1);
$count1= scalar(grep {defined $_} @file_list1), "\n"; 
for(my $i=0;$i<=$count1;$i++)
{

	if($file_list1[$i] =~ m/^[0-9a-zA-Z].*/)
	{
		push(@files1, $file_list1[$i]);
	}
}

open(FILE1,<@ARGV[0]>);
@array1=<FILE1>;

$j1=0;
foreach my $k(@array1)
{
	if($k=~/ : /)
	{
		chop($k);
		@a=split(/ : /,$k);
		$b[$j1]=$a[1];
		$j1++;
	}
}	
close FILE1;

$a1=1;
foreach(@files1)
{	
	print "file:$a1  ";
	system("fastx_trimmer -Q33 $b[0] -i ./raw_data/$_ | fastx_trimmer -Q33 $b[1]| fastq_quality_trimmer -Q33 $b[2] $b[3]|fastq_quality_filter -Q33 $b[4] $b[5] -o ./cleaned_data/$_");
	print "done  ";
	$a1++;
}
system("mkdir concatenated_files");
foreach $cleaned_file (@files1)
{
	@temp_array1 = split(/_R/,$cleaned_file);
	if($temp_array1[1] =~ m/1_001.fastq/)
	{
		@temp_array2=($temp_array1[0],'_R2_001.fastq'); 
		my $read2_file=join("",@temp_array2);
		system("cat  ./cleaned_data/$cleaned_file  ./cleaned_data/$read2_file >> ./concatenated_files/$temp_array1[0].fastq");
	}
}
system("rm -r cleaned_data");
###############################
opendir(DIR2,"genes");
@genes_files=readdir(DIR2);
closedir(DIR2);

$count_genes=scalar(grep {defined $_} @genes_files), "\n";
for ($genes=0;$genes<=$count_genes;$genes++)
{
	if( $genes_files[$genes] =~ m/^[0-9a-zA-Z].*/ )
	{
		push(@genes,$genes_files[$genes]);
	}
}
foreach my $gene (@genes)
{
	
               chdir "genes" or die;
	       open(FILE2,"$gene");
               my $c= <FILE2>;
               $/=">";
               my $d= <FILE2>;
               my (@genebp) = split /\n/, $d;
               my $genelength = join('',@genebp);
	       $gene_num=length($genelength);
	       push(@names,$gene,$gene_num);	
               $/="\n";
               close FILE2;
	       chdir "./../" or die;
}

%hash=@names;
opendir(DIR3,"concatenated_files");
@samples_files=readdir(DIR3);
closedir(DIR3);
$count_samples=scalar(grep {defined $_} @samples_files), "\n";
for($samples=0;$samples<=$count_samples;$samples++)
{
	if( $samples_files[$samples] =~ m/^[0-9a-zA-Z].*/ )
	{	
		push(@samples,$samples_files[$samples]);
	}
}

system("mkdir ./consensus");

	open(FILE3,Makefile);
	@array2=<FILE3>;
	close FILE3;

$sample_number=1;
foreach my $sample (@samples)
{
	system ("mkdir sample_$sample");
	chdir "./sample_$sample" or die;
	@name_split=split(/\./,$sample);
	my $filename="$name_split[0]";
	open(FILE5,">$filename");
	print color 'bold blue';	
	print "Analysing sample_number:$sample_number\n\n";
	print color 'reset';

	foreach my $gene (@genes)
	{
	  	system("mkdir gene_$gene");
		chdir "./gene_$gene" or die;		
		
		foreach my $i (@array2)
		{
			if($i =~ m/READS=/)
			{
				my $a="./../../concatenated_files/$sample";
				$i = "READS=$a";
			}
			if($i =~ m/TEMPLATE=/)
			{
				my $b="./../../genes/$gene";
				$i = "TEMPLATE=$b";
			}
			
				open(FILE6,'>>Makefile');
				chomp($i);
				print FILE6 "$i\n";
				close FILE6;
		}
	
	system ("make single_step TYPE=solexa ORIENT=linear PID=same");
		
	open(FILE4,'Final_Assembly');
	my $juncy=<FILE4>;
	$/ = ">";
	my $junk = <FILE4>;
	chomp($junk);
	my (@seqlines) = split /\n/, $junk;
	my $sequence = join('',@seqlines);
	my $seq_length= length($sequence);
	$/="\n";

	if( ($seq_length/$hash{$gene}) >= $length_criteria )
	{  
		print color 'magenta';
		print "\nConsensus created with a 95% identity to referance!\n\n";
		print color 'reset';

		@list = split(/\./, $gene);
		print FILE5 ">$list[0]\n";		
        	print FILE5 "$junk\n";
	}
	else
	{
		system("make clean");
		system("make single_step TYPE=solexa ORIENT=linear PID=medium");

		open(FILE4,'Final_Assembly');
        	my $juncy=<FILE4>;
        	$/ = ">";
        	my $junk = <FILE4>;
        	chomp($junk);
		my (@seqlines) = split /\n/, $junk;
        	my $sequence = join('',@seqlines);
        	my $seq_length= length($sequence);
		$/="\n";		

		if( ($seq_length/$hash{$gene}) >= $length_criteria )
        	{
			print color 'magenta';
                	print "\nConsensus created with an 85% identity to referance!\n\n";
			print color 'reset';

			@list= split (/\./,$gene);
                	print FILE5 ">$list[0]\n";
        		print FILE5 "$junk\n";
		}
		else
		{
			system("make clean");
			system("make single_step TYPE=solexa ORIENT=linear PID=desperate");	
			open(FILE4,'Final_Assembly');
                	my $juncy=<FILE4>;
                	$/ = ">";
                	my $junk = <FILE4>;
                	chomp($junk);
                	my (@seqlines) = split /\n/, $junk;
                	my $sequence = join('',@seqlines);
                	my $seq_length= length($sequence);
			$/="\n";

			if( ($seq_length/$hash{$gene}) >= $length_criteria )
			{
				print color 'magenta';
				print "\nConsensus created with less than 85% identity to referance!\n\n";
				print color 'reset';

				($s1, $s2) = (split/\./,$gene);
                		print FILE5 ">$s1\n";
				print FILE5 "$junk\n";
			}
			else
			{
				system("make clean");
			}
		}			
	}
	chdir "./../" or die;
	}
	$sample_number++;
	system("mv $filename ./../consensus");		
	chdir "./../" or die;
	system("rm -r sample_*");
}
close FILE4;
close FILE5;
system("rm -r concatenated_files");
###############################################
opendir(DIR4,"consensus");
@consensus_files=readdir(DIR4);
closedir(DIR4);
$count2= scalar(grep {defined $_} @consensus_files); 
system("mkdir output");
open(HEADER,'<header');
while(<HEADER>)
{
	my $line=$_;
	chop($line);
	@header_array=split(/\t/,"$line");
}
close HEADER;

for(my $i=0;$i<=$count2;$i++)
{

	if($consensus_files[$i] =~ m/^[0-9a-zA-Z].*/)
	{
		push(@consensus_file, $consensus_files[$i]);
	}
}
system("mkdir paralog_warnings");
chdir "./consensus" or die;

foreach my $file (@consensus_file)
{                                    
	print color 'magenta';
	print "Sample consensus file:$file\n";
	print color 'reset';
	open(READ,"$file") || die ("Error opening $file $!");
	open(WRITE,">../output/$file") || die ("Error opening $file $!");
	
	while (my $line = <READ>)
	{
		chomp($line);
		if ($line=~m/^>/) 
		{ 
			print WRITE "\n",$line,"\n"; 
		}
		else 
		{ 	
			print WRITE "$line"; 
		}
	}
	print WRITE "\n";
	close READ; 
		
	open(WRITE,"../output/$file") || die ("Error opening $file $!");
	while(my $line = <WRITE>)
	{
		chomp($line);	
		if($line =~ m/^(A|T|G|C)/ )	
                {
			@linelength=split(//,$line);
			$seqlength=@linelength;
			my $seq = $line;
			my $url="http://127.0.0.1/cgi-bin/bigsdb/bigsdb.pl?locus=0&order=locus&sequence=".$seq."&accession=&submit=Submit&db=bigsdb_listeria_seqdef&page=sequenceQuery";
			my $content =get $url;
			open(OUT,'>>outfile');
			print OUT "$content\n";
			
			system("cat outfile|wc -l >> numfile"); 
			system ("rm numfile");                   
			open(GTOUT,'<outfile');     
			@newarray=<GTOUT>;
			$totaline=@newarray; 
			for($j=0;$j<=$totaline;$j++)
			{
				if($newarray[$j] =~ m/^1 exact match found*/)
				{
					print "1 exact match found!\n";
					system("grep '<table' outfile > out");
					open(OUT1,out);
					@tableline=<OUT1>;
					chop($tableline[0]);
					@firstarray=split(/=/,$tableline[0]);
					@secondarray=split(/ \(/,$firstarray[7]);
					@thirdarray=split(/">/,$secondarray[0]);
					@fourarray=split(/ /,$thirdarray[1]);
					if( $fourarray[1] =~ m/^[A-Za-z0-9]/ )
					{
						 $loci_name=join("_","$fourarray[0]","$fourarray[1]");
					}
					elsif( $fourarray[1] == "")
					{
						$loci_name=$fourarray[0];
					}
					print "$loci_name::$thirdarray[0]\n";
					system("rm out");
				}
				elsif($newarray[$j] =~ m/^2 exact matches found*/)
				{
					system("grep 'td1' outfile > first");
					system("grep 'td2' outfile > second");
					open(OUT1,first);
					@firstline=<OUT1>;
					chop($firstline[0]);
					open(OUT2,second);
					@secondline=<OUT2>;
					chop($secondline[0]);
					@firsthit1=split(/=/,$firstline[0]);
					@secondhit1=split(/=/,$secondline[0]);

					@check1=split(/<\/td><td>/,$firsthit1[7]);
                			@check2=split(/<\/td><td>/,$secondhit1[6]);

					if ($check1[1] ==$check2[1] && $check1[2]==$check2[2] && $check1[2]==1 && $check2[1]==$seqlength)
	   				{		
						print "2 exact matches found!\n";
						@firsthit2=split(/ \(/,$firsthit1[7]);
						@firsthit3=split(/">/,$firsthit2[0]);
						@firsthit4=split(/ /,$firsthit3[1]);
						if( $firsthit4[1] =~ m/^[A-Za-z0-9]/ )
						{
							$lociname1=join("_","$firsthit4[0]","$firsthit4[1]");
						}
						elsif( $firsthit4[1] == "")
						{
							$lociname1=$firsthit4[0];
						}
											
						@secondhit2=split(/">/,$secondhit1[6]);
						@secondhit3=split(/ \(/,$secondhit2[1]);
						@secondhit4=split(/ /,$secondhit3[0]);
						if( $secondhit4[1] =~ m/^[A-Za-z0-9]/ )
						{
							$lociname2=join("_","$secondhit4[0]","$secondhit4[1]");
							$loci_name=0;
						}
						elsif( $secondhit4[1] == "")
						{
							$lociname2=$secondhit4[0];
							$loci_name=0;
						}
						open(WARNING,">>../paralog_warnings/$file");
						print WARNING "$lociname1:$firsthit3[0] matches $lociname2:$secondhit2[0]\n";
						print "$lociname1:$firsthit3[0] matches $lociname2:$secondhit2[0]\n";
					}
					elsif($check1[2]==1 && $check1[1]==$seqlength)
	   				{   
						print "2 exact matches found but only the 1st is a perfect match\n";
						@firsthit2=split(/ \(/,$firsthit1[7]);
	    	           		        @firsthit3=split(/">/,$firsthit2[0]);
        	        			@firsthit4=split(/ /,$firsthit3[1]);
                				if( $firsthit4[1] =~ m/^[A-Za-z0-9]/ )
                				{
                        				$lociname1=join("_","$firsthit4[0]","$firsthit4[1]");
                				}
                				elsif( $firsthit4[1] == "")
                				{
                        				$lociname1=$firsthit4[0];
                				}
						print "$lociname1::$firsthit3[0]\n";
	    				}
	    				elsif($check2[2]==1 && $check2[2]==$seqlength)
	    				{ 
						print "2 exact matches found but only the 2nd is a perfect match\n";
						@secondhit2=split(/">/,$secondhit1[6]);
                				@secondhit3=split(/ \(/,$secondhit2[1]);
                				@secondhit4=split(/ /,$secondhit3[0]);
                				if( $secondhit4[1] =~ m/^[A-Za-z0-9]/ )
                				{
                        				$lociname2=join("_","$secondhit4[0]","$secondhit4[1]");
                       	        			$loci_name=0;
                				}
                				elsif( $secondhit4[1] == "")
                				{
                        				$lociname2=$secondhit4[0];
                        				$loci_name=0;
                				}
						print "$lociname2::$secondhit2[0]\n";
	     				}
					system("rm first second");
				}
				elsif($newarray[$j] =~ m/^<p style="margin-top:0.5em">Closest match:/)
				{
					print "closest match found!\n";
               			        @tableline=split(/:/,$newarray[$j+2]);
			                @firstarray=split(/</,$tableline[1]);
                			$thirdarray[0] = join(//,"closest match:","$firstarray[0]");
                			@secondarray=split(/ \(/,$tableline[0]);
                			@fourarray=split(/ /,$secondarray[0]);
                			if( $fourarray[1] != "" )
                			{
                         			$loci_name=join("_","$fourarray[0]","$fourarray[1]");
                			}
                			elsif( $fourarray[1] == "")
                			{
                        			$loci_name=$fourarray[0];
                			}
					 print "$loci_name::$thirdarray[0]\n";
				}
			}
			system("rm outfile");
					
                        close GTOUT;
			my $c=0;
			foreach my $value (@header_array)
			{
				if($value eq $loci_name)
				{
					$cell[$c]="$thirdarray[0]";
				}
				elsif($value eq $lociname1)
				{
					$cell[$c]="$firsthit3[0]";
				}
				elsif($value eq $lociname2)
				{
					$cell[$c]="$secondhit2[0]";
				}
			    	$c++;
			}
		}
		@tableline = (); @firstline=(); @secondline=(); 
		@firstarray = (); @firsthit1=(); @firsthit2=(); @firsthit3=(); @firsthit4=();
		@secondarray = (); @secondhit1=(); @secondhit2=(); @secondhit3=(); @secondhit4=();
		@thirdarray = (); @fourarray=();
		$loci_name =""; 
		$lociname1=""; 
		$lociname2="";
	}
	 my $sample = join("\t",@cell);
	 sleep(5);
	 @cell=();
	 @name_array = split(/\./,$file); 	
	 my $final_sample = "$name_array[0] $sample";
	 close WRITE;
         open(HEADER,'>>../header');
         print HEADER "$final_sample\n";
	 sleep(2);
         close HEADER;
	 system("rm  ../output/$file");
	 sleep(20);
}
system("rm -r ../output");
system("cd ./../");

