=begin JA

= ����ӥ�إ��ۥ���԰���μ¸� 

# * �����̰�ϯ (sugiyama)
#   * $Id: catseye.rd,v 1.1 2014/02/26 07:09:53 sugiyama Exp $

=end JA
=begin EN

= KelvinHelmholtz Instability Experiment 

# * Ko-ichiro Sugiyama
#   * $Id: catseye.rd,v 1.1 2014/02/26 07:09:53 sugiyama Exp $

=end EN

=begin JA
����ӥ�إ��ۥ���԰���μ¸���¹Ԥ�����ˡ���������ޤ�. 
=end JA
=begin EN
A method to perform a KelvinHelmholtz instability experiment is described. 
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
�ʲ��Ǥ� deepconv/arare5 �������ǥ��쥯�ȥ���٤� ../deepconv-exp/catseye-exp �ǥ��쥯�ȥ�������, �����Ǽ¸���Ԥ����Ȥˤ��ޤ�. 
���Τ褦��  ../deepconv-exp/catseye-exp �ǥ��쥯�ȥ�������, ��ư���Ƥ�������. 

  $ mkdir -p ../deepconv-exp/catseye-exp
  $ cd ../deepconv-exp/catseye-exp

����, ���Υǥ��쥯�ȥ�˼¹ԥե����������ե������֤����������ޤ�. 

  $ mkdir bin
  $ mkdir conf

�Ǹ�˺��������ǥ��쥯�ȥ�˼¹ԥե����������ե�����򥳥ԡ����ޤ�. 

  $ cp ../../arare5-YYYYMMDD/src/main/arare bin
  $ cp ../../arare5-YYYYMMDD/src/main/arare_init-data bin
  $ cp ../../arare5-YYYYMMDD/exp_setup_files/*catseye*.conf conf

�ʤ�, �¹ԥե����������ե����� (NAMELIST �ե�����) �������, �ɤΥǥ��쥯�ȥ�ˤ����Ƥ�׻���Ԥ����Ȥ��Ǥ��ޤ�. 


=end JA
=begin EN

== Preparation of a directory for experiments

Let us move to the top directory of deepconv/arare5 src tree
(assuming arare5-YYYYMMDD in the following example). 
Here, we perform an experiment in ../arare5-exp/catseye-exp directory. 
Please create the directory and enter there as follows: 

  $ mkdir -p ../deepconv-exp/catseye-exp
  $ cd ../deepconv-exp/catseye-exp

Then, pleaase create the directories for executable files and configuration
files as follows:

  $ mkdir bin
  $ mkdir conf

Finally, executable files and configuration files are copied as follows: 

  $ cp ../../arare5-YYYYMMDD/src/main/arare bin
  $ cp ../../arare5-YYYYMMDD/src/main/arare_init-data bin
  $ cp ../../arare5-YYYYMMDD/exp_setup_files/*catseye*.conf conf

Note that you can perform an experiment in any directory by using executable files and configuration (NAMELIST) files.

=end EN


=begin JA

== ����ͥǡ����ե�����κ���

arare_init-data �� 
((<arare-catseye_init-data.conf|URL:../../exp_setup_files/arare-catseye_init-data.conf>))
���Ѥ��ƽ���ͥե����� arare-catseye_init.nc ��������ޤ�. 

  $ ./bin/arare_init-data -N=conf/arare-catseye_init-data.conf

   *** MESSAGE [main] ***  Namelist file is 'arare-catseye_init-data.conf'
   *** MESSAGE [argset_init] ***  NAMELIST FILE = arare-catseye_init-data.conf
   *** MESSAGE [namelist_util] ***  ----- Initialization Messages -----
   *** MESSAGE [namelist_util] ***    MaxNmlArySize = 256
   *** MESSAGE [namelist_util] ***  -- version = $Name:  $$Id: catseye.rd,v 1.1 2014/02/26 07:09:53 sugiyama Exp $
   *** MESSAGE [gridset_init] ***  xsub = 1
                         : 
   *** MESSAGE [main] ***  Making Initial data....
   *** MESSAGE [main] ***  Making Initial data (basic) named IsoThermal...
   *** MESSAGE [main] ***  Making Initial data (disturb) named tanh...
   *** MESSAGE [main] ***  Making Initial wind data (disturb) named tanh...
   *** MESSAGE [main] ***  Output variables into netCDF file...
   *** MESSAGE [restartfileioIO_init] ***  InputFile  =
   *** MESSAGE [restartfileioIO_init] ***  OutputFile = arare-catseys_init.nc


== �¸��μ¹�

�¹ԥե����� arare �� NAMELIST �ե�����
((<arare-catseye.conf|URL:../../exp_setup_files/arare-catseye.conf>))
���Ѥ���, �ʲ��Τ褦�� arare ��¹Ԥ��Ƥ�������. 
�ץ���ཪλ�ˤϿ�ʬ������ޤ�. 

(�ʤ�, ��������ѥ���Ķ��Ǥϰʲ�����ˡ�ǥץ�����
�¹Ԥ��뤳�ȤϤǤ��ʤ��Τ���դ��Ƥ�������. ���ξ��μ¹���ˡ
�˴ؤ��Ƥ�, ���δĶ��ǤΥץ����¹ԥޥ˥奢��ʤɤ򻲾Ȥ�������. )

  $ ./bin/arare -N=conf/arare-catseye.conf | tee catseye.log

   *** MESSAGE [main] ***  Namelist file is 'arare-catseye.conf'
   *** MESSAGE [argset_init] ***  NAMELIST FILE = arare-catseye.conf
   *** MESSAGE [namelist_util] ***  ----- Initialization Messages -----
   *** MESSAGE [namelist_util] ***    MaxNmlArySize = 256
   *** MESSAGE [namelist_util] ***  -- version = $Name:  $$Id: catseye.rd,v 1.1 2014/02/26 07:09:53 sugiyama Exp $
   *** MESSAGE [timeset_init] ***  DelTimeLong  = 1.
                     :
   *** MESSAGE [HistoryClose] ***  "arare-catseye_Exner.nc" is closed
   *** MESSAGE [HistoryClose] ***  "arare-catseye_Km.nc" is closed
    
   ############## CPU TIME SUMMARY ################
   initialization         0.320010E-01
   time-integration       0.610998E+02  (1.02 minutes)
   ------------------------------------------------
          TOTAL TIME =    0.611318E+02  (1.02 minutes)

���ξ��, �� 1 ʬ�λ�����ʬ���פ��ޤ���. 
�׻���̤� VelX.nc �� PTemp.nc �Ȥ��ƽ��Ϥ���ޤ�. 
�ޤ�, �ꥹ�����ȥե����뤬 arare-catseye_restart.nc �Ȥ��ƽ��Ϥ���ޤ�. 

=end JA

=begin EN

== Create initial data file

Create initial data file "denscurrent-dry_restart.nc"
using "bin/arare_init-data" and ((<"arare-catseye_init-data.conf"|URL:../../exp_setup_files/arare-catseye_init-data.conf>)).
  
  $ ./bin/arare_init-data -N=conf/arare-catseye_init-data.conf

   *** MESSAGE [main] ***  Namelist file is 'arare-catseye_init-data.conf'
   *** MESSAGE [argset_init] ***  NAMELIST FILE = arare-catseye_init-data.conf
   *** MESSAGE [namelist_util] ***  ----- Initialization Messages -----
   *** MESSAGE [namelist_util] ***    MaxNmlArySize = 256
   *** MESSAGE [namelist_util] ***  -- version = $Name:  $$Id: catseye.rd,v 1.1 2014/02/26 07:09:53 sugiyama Exp $
   *** MESSAGE [gridset_init] ***  xsub = 1
                         : 
   *** MESSAGE [main] ***  Making Initial data....
   *** MESSAGE [main] ***  Making Initial data (basic) named IsoThermal...
   *** MESSAGE [main] ***  Making Initial data (disturb) named tanh...
   *** MESSAGE [main] ***  Making Initial wind data (disturb) named tanh...
   *** MESSAGE [main] ***  Output variables into netCDF file...
   *** MESSAGE [restartfileioIO_init] ***  InputFile  =
   *** MESSAGE [restartfileioIO_init] ***  OutputFile = arare-catseys_init.nc

== Run the experiment

Using an executable files 'arare' and a NAMELIST file
((<arare-catseye.conf|URL:../../exp_setup_files/arare-catseye.conf>)), 
execute 'arare' as follows. 
This program will be finished in few minutes - tens of minutes. 

  $ ./bin/arare -N=conf/arare-catseye.conf | tee catseye.log

   *** MESSAGE [main] ***  Namelist file is 'arare-catseye.conf'
   *** MESSAGE [argset_init] ***  NAMELIST FILE = arare-catseye.conf
   *** MESSAGE [namelist_util] ***  ----- Initialization Messages -----
   *** MESSAGE [namelist_util] ***    MaxNmlArySize = 256
   *** MESSAGE [namelist_util] ***  -- version = $Name:  $$Id: catseye.rd,v 1.1 2014/02/26 07:09:53 sugiyama Exp $
   *** MESSAGE [timeset_init] ***  DelTimeLong  = 1.
                     :
   *** MESSAGE [HistoryClose] ***  "arare-catseye_Exner.nc" is closed
   *** MESSAGE [HistoryClose] ***  "arare-catseye_Km.nc" is closed
    
   ############## CPU TIME SUMMARY ################
   initialization         0.320010E-01
   time-integration       0.610998E+02  (1.02 minutes)
   ------------------------------------------------
          TOTAL TIME =    0.611318E+02  (1.02 minutes)

In this case, about 1 minites integration is performed. 
History data are output to 'VelX.nc' and 'PTemp.nc' etc., 
and a restart data is output to 'arare-catseye_restart.nc'. 

=end EN

=begin JA
== ��̤βĻ벽
((<��ñ�ʲ��ϡ��Ļ벽|URL:./visualization.htm>)) �򻲾Ȥ��Ƥ�������. 

((<"IMG:catseye_img01.png">))

=end JA
=begin EN
== Visualization

Please see ((<First step analysis and visualization|URL:./visualization.htm.en>)).

((<"IMG:catseye_img01.png">))

=end EN


=begin JA
== ����ʸ��
* ���ڤȺ縶, CReSS �ޥ˥奢���� 2 ��.
=end JA

=begin EN
== References
* Tsuboki and Sakakibara, CReSS User Manual (in Japanese)
=end EN




=begin HTML
<hr />
<small>
  $Id: catseye.rd,v 1.1 2014/02/26 07:09:53 sugiyama Exp $
</small>
=end HTML

