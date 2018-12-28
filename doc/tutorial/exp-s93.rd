=begin JA

= ����ή�μ¸� by Straka et al. (1993)

# * �����̰�ϯ (sugiyama)
#   * $Id: exp-s93.rd,v 1.2 2014/03/04 08:09:19 sugiyama Exp $

=end JA
=begin EN

= Density Current Experiment by Straka et al. (1993)

# * Yoshiyuki O. Takahashi (yot), Shin-ichi Takehiro (takepiro)
#   * $Id: exp-s93.rd,v 1.2 2014/03/04 08:09:19 sugiyama Exp $

=end EN

=begin JA
Straka et al. (1993) �ǹԤ�줿, ����ή�μ¸���¹Ԥ�����ˡ���������ޤ�. 
=end JA
=begin EN
A method to perform a density current experiment by Straka et al. (1993) is described. 
=end EN


=begin JA
== ����
�ܼ¸��ϰʲ��� 3 �ĤΥ��ƥåפǹԤ��ޤ�.

  * �¸��ǥ��쥯�ȥ�ν���
  * ����ͤν���
  * �¸��μ¹�
=end JA
=begin EN
== Overview
This experiment is performed with the following 3 steps:

  * Preparation of directory for experiments
  * Preparation of initial condition
  * Execution of experiments
=end EN


=begin JA

== �¸��ѥǥ��쥯�ȥ����

deepconv �ο��ͼ¸��ϥ������ĥ꡼�����ǤϹԤ鷺
�������ĥ꡼�Ȥ��̤γ����ǥ��쥯�ȥ�ˤƹԤ����Ȥ�侩�������ޤ�. 

�ޤ� deepconv/arare5 �������Υȥåץǥ��쥯�ȥ�(�ʲ�����Ǥ� arare5-YYYYMMDD �Ȥ���)�˰�ư���Ƥ�������. 
�ʲ��Ǥ� deepconv/arare5 �������ǥ��쥯�ȥ���٤� ../deepconv-exp/s93-exp �ǥ��쥯�ȥ�������, �����Ǽ¸���Ԥ����Ȥˤ��ޤ�. 
���Τ褦��  ../deepconv-exp/s93-exp �ǥ��쥯�ȥ�������, ��ư���Ƥ�������. 

  $ mkdir -p ../deepconv-exp/s93-exp
  $ cd ../deepconv-exp/s93-exp

����, ���Υǥ��쥯�ȥ�˼¹ԥե����������ե������֤����������ޤ�. 

  $ mkdir bin
  $ mkdir conf

�Ǹ�˺��������ǥ��쥯�ȥ�˼¹ԥե����������ե�����򥳥ԡ����ޤ�. 

  $ cp ../../arare5-YYYYMMDD/src/main/arare bin
  $ cp ../../arare5-YYYYMMDD/src/main/arare_init-data bin/
  $ cp ../../arare5-YYYYMMDD/exp_setup_files/arare-DensCurrent-dry*.conf conf/

�ʤ�, �¹ԥե����������ե����� (NAMELIST �ե�����) �������, �ɤΥǥ��쥯�ȥ�ˤ����Ƥ�׻���Ԥ����Ȥ��Ǥ��ޤ�. 


=end JA
=begin EN

== Preparation of a directory for experiments

Let us move to the top directory of deepconv/arare5 src tree
(assuming arare5-YYYYMMDD in the following example). 
Here, we perform an experiment in ../arare5-exp/s93-exp directory. 
Please create the directory and enter there as follows: 

  $ mkdir -p ../deepconv-exp/s93-exp
  $ cd ../deepconv-exp/s93-exp

Then, pleaase create the directories for executable files and configuration
files as follows:

  $ mkdir bin
  $ mkdir conf

Finally, executable files and configuration files are copied as follows: 

  $ cp ../../arare5-YYYYMMDD/src/main/arare bin
  $ cp ../../arare5-YYYYMMDD/src/main/arare_init-data bin/
  $ cp ../../arare5-YYYYMMDD/exp_setup_files/arare-DensCurrent-dry*.conf conf/

Note that you can perform an experiment in any directory by using executable files and configuration (NAMELIST) files.

=end EN


=begin JA

== ����ͥǡ����ե�����κ���

arare_init-data �� 
((<arare-DensCurrent-dry_init.conf|URL:../../exp_setup_files/arare-DensCurrent-dry_init.conf>))
���Ѥ��ƽ���ͥե����� denscurrent-dry_restart.nc ��������ޤ�. 

  $ ./bin/arare_init-data -N=conf/arare-DensCurrent-dry_init.conf

   *** MESSAGE [main] ***  Namelist file is 'arare-DensCurrent-dry_init.conf'
   *** MESSAGE [argset_init] ***  NAMELIST FILE = arare-DensCurrent-dry_init.conf
   *** MESSAGE [namelist_util] ***  ----- Initialization Messages -----
   *** MESSAGE [namelist_util] ***    MaxNmlArySize = 256
   *** MESSAGE [namelist_util] ***  -- version = $Name:  $$Id: exp-s93.rd,v 1.2 2014/03/04 08:09:19 sugiyama Exp $
   *** MESSAGE [gridset_init] ***  xsub = 1
                         : 
   *** MESSAGE [main] ***  Making Initial data....
   *** MESSAGE [main] ***  Making Initial data (basic) named DRY...
   *** MESSAGE [main] ***  Making Initial data (disturb) named CosXZ...
   *** MESSAGE [main] ***  Output variables into netCDF file...
   *** MESSAGE [restartfileioIO_init] ***  InputFile  =
   *** MESSAGE [restartfileioIO_init] ***  OutputFile = denscurrent-dry_restart.nc

== �¸��μ¹�

�¹ԥե����� arare �� NAMELIST �ե�����
((<arare-DensCurrent-dry.conf|URL:../../exp_setup_files/arare-DensCurrent-dry.conf>))
���Ѥ���, �ʲ��Τ褦�� arare ��¹Ԥ��Ƥ�������. 
�ץ���ཪλ�ˤϿ�ʬ������ʬ������ޤ�. 

(�ʤ�, ��������ѥ���Ķ��Ǥϰʲ�����ˡ�ǥץ�����
�¹Ԥ��뤳�ȤϤǤ��ʤ��Τ���դ��Ƥ�������. ���ξ��μ¹���ˡ
�˴ؤ��Ƥ�, ���δĶ��ǤΥץ����¹ԥޥ˥奢��ʤɤ򻲾Ȥ�������. )

  $ ./bin/arare -N=conf/arare-DensCurrent-dry.conf | tee s93.log

   *** MESSAGE [main] ***  Namelist file is 'arare-DensCurrent-dry.conf'
   *** MESSAGE [argset_init] ***  NAMELIST FILE = arare-DensCurrent-dry.conf
   *** MESSAGE [namelist_util] ***  ----- Initialization Messages -----
   *** MESSAGE [namelist_util] ***    MaxNmlArySize = 256
   *** MESSAGE [namelist_util] ***  -- version = $Name:  $$Id: exp-s93.rd,v 1.2 2014/03/04 08:09:19 sugiyama Exp $
   *** MESSAGE [timeset_init] ***  DelTimeLong  = 1.
                     :
   *** MESSAGE [HistoryClose] ***  "denscurrent-dry_ElstEnrgy.nc" is closed
   *** MESSAGE [HistoryClose] ***  "denscurrent-dry_PotEnrgy.nc" is closed
    
   ############## CPU TIME SUMMARY ################
   initialization         0.720040E-01
   time-integration       0.234947E+03  (3.92 minutes)
   ------------------------------------------------
          TOTAL TIME =    0.235019E+03  (3.92 minutes)

���ξ��, �� 4 ʬ�λ�����ʬ���Ԥ��ޤ�. 
�׻���̤� VelX.nc �� PTemp.nc �Ȥ��ƽ��Ϥ���ޤ�. 
�ޤ�, �ꥹ�����ȥե����뤬 denscurrent-dry_restart2.nc �Ȥ��ƽ��Ϥ���ޤ�. 

=end JA

=begin EN

== Create initial data file

Create initial data file "denscurrent-dry_restart.nc"
using "bin/arare_init-data" and ((<"arare-DensCurrent-dry_init.conf"|URL:../../exp_setup_files/arare-DensCurrent-dry_init.conf>)).
  
  $ ./bin/arare_init-data -N=conf/arare-DensCurrent-dry_init.conf

   *** MESSAGE [main] ***  Namelist file is 'arare-DensCurrent-dry_init.conf'
   *** MESSAGE [argset_init] ***  NAMELIST FILE = arare-DensCurrent-dry_init.conf
   *** MESSAGE [namelist_util] ***  ----- Initialization Messages -----
   *** MESSAGE [namelist_util] ***    MaxNmlArySize = 256
   *** MESSAGE [namelist_util] ***  -- version = $Name:  $$Id: exp-s93.rd,v 1.2 2014/03/04 08:09:19 sugiyama Exp $
   *** MESSAGE [gridset_init] ***  xsub = 1
                         : 
   *** MESSAGE [main] ***  Making Initial data....
   *** MESSAGE [main] ***  Making Initial data (basic) named DRY...
   *** MESSAGE [main] ***  Making Initial data (disturb) named CosXZ...
   *** MESSAGE [main] ***  Output variables into netCDF file...
   *** MESSAGE [restartfileioIO_init] ***  InputFile  =
   *** MESSAGE [restartfileioIO_init] ***  OutputFile = denscurrent-dry_restart.nc

== Run the experiment

Using an executable files 'arare' and a NAMELIST file
((<arare-DensCurrent-dry.conf|URL:../../exp_setup_files/arare-DensCurrent-dry.conf>)), 
execute 'arare' as follows. 
This program will be finished in few minutes - tens of minutes. 

  $ ./bin/arare -N=conf/arare-DensCurrent-dry.conf | tee s93.log

   *** MESSAGE [main] ***  Namelist file is 'arare-DensCurrent-dry.conf'
   *** MESSAGE [argset_init] ***  NAMELIST FILE = arare-DensCurrent-dry.conf
   *** MESSAGE [namelist_util] ***  ----- Initialization Messages -----
   *** MESSAGE [namelist_util] ***    MaxNmlArySize = 256
   *** MESSAGE [namelist_util] ***  -- version = $Name:  $$Id: exp-s93.rd,v 1.2 2014/03/04 08:09:19 sugiyama Exp $
   *** MESSAGE [timeset_init] ***  DelTimeLong  = 1.
                     :
   *** MESSAGE [HistoryClose] ***  "denscurrent-dry_ElstEnrgy.nc" is closed
   *** MESSAGE [HistoryClose] ***  "denscurrent-dry_PotEnrgy.nc" is closed
    
   ############## CPU TIME SUMMARY ################
   initialization         0.720040E-01
   time-integration       0.234947E+03  (3.92 minutes)
   ------------------------------------------------
          TOTAL TIME =    0.235019E+03  (3.92 minutes)

In this case, about 4 minites integration is performed. 
History data are output to 'VelX.nc' and 'PTemp.nc' etc., 
and a restart data is output to 'denscurrent-dry_restart2.nc'. 

=end EN

=begin JA
== ��̤βĻ벽
((<��ñ�ʲ��ϡ��Ļ벽|URL:./visualization.htm>)) �򻲾Ȥ��Ƥ�������. 

((<"IMG:exp-s93_img01.png">))

=end JA
=begin EN
== Visualization

Please see ((<First step analysis and visualization|URL:./visualization.htm.en>)).

((<"IMG:exp-s93_img01.png">))

=end EN


=begin JA
== ����ʸ��
=end JA

=begin EN
== References
=end EN

=begin
* Straka, J. M., Wilhelmson, R. B., Wicker, A. L., Anderson, J. R., Droegemeier, K. K., 1993 : Numerical solutions of a non-linear density current : A benchmark solution and comparisons. International Journal for Numerical Methods in Fluids , 17, 1--22

=end


=begin HTML
<hr />
<small>
  $Id: exp-s93.rd,v 1.2 2014/03/04 08:09:19 sugiyama Exp $
</small>
=end HTML

