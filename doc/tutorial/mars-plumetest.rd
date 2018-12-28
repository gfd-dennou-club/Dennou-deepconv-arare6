=begin TOPLINK
[((<English|URL:mars-plumetest.htm.en>)) |
((<Japanese|URL:mars-plumetest.htm>))]
[((<GFD Dennou Club|URL:http://www.gfd-dennou.org>)) |
((<deepconv Project|URL:http://www.gfd-dennou.org/library/deepconv>))]
=end TOPLINK


=begin JA

= ���������ޥ�徺�¸�

#* ����
#  * 2011/03/01(����ã��) ��������


# * ���� ã�� (yamasita), ���� �̰�ϯ (sugiyama), ���� ���� (odakker)
#   * $Id: mars-plumetest.rd,v 1.1 2011/06/21 15:42:02 sugiyama Exp $

=end JA
=begin EN

= A Numerical experiment of ascending hot plume

#* History
#  * 2011/03/01 (Tatsuya Yamashita) Initial release


# * Tatsuya YAMASHITA (yamasita), Ko-ichiro SUGIYAMA (sugiyama), Masatsugu ODAKA (odakker)
#   * $Id: mars-plumetest.rd,v 1.1 2011/06/21 15:42:02 sugiyama Exp $

=end EN


=begin JA

�����Ǥ�, ����¾(2006)�ǹԤʤ�줿, 
ñ�ȤΥۥåȥץ�塼���徺������¸���¹Ԥ�����ˡ���������ޤ�. 

���η׻��ˤϰʲ���ʪ���������Ѥ��Ƥ��ޤ�. 

  * ��ή�Ȼ�
    * 1.5 ���Υ������㡼��ǥ� (Klemp and Wilhelmson, 1978)
  * �ŷ롦��ȯ
    * �Ȼ���Ĺ (e.g., Tobie et al., 2003)

=end JA
=begin EN

A method to perform an test experiment
of an ascending hot plume is described. 

Following physical processes are used in this experiment.

  * Turbulent mixing
    * 1.5 order closure (Klemp and Wilhelmson, 1978)
  * Condensation and evaporation
    * Diffusional growth (e.g., Tobie et al., 2003)

=end EN


=begin JA
== ����

�ܼ¸��ϥǥե���Ȥ�����ǹԤʤ����Ȥ�����ޤ�. 
�ǥե��������Ǥη׻��򤹤���Ϥ��Ρֳ��ספΤߤ��ɤ�н�ʬ�Ǥ�. 
�ǥե���ȤȤϰۤʤ�����Ƿ׻��򤹤����, �־ܺ١פ�����������. 


��������ʬ�ŷ���ή�ѤΥ������Υ���ѥ��뤬��λ�����, 
�׻���¹Ԥ��뤳�Ȥ��Ǥ��ޤ�. 
����ѥ�����ˤ� make �ǤϤʤ�, make mmc ���ǤĤ��Ȥ���դ��Ʋ�����. 
�㤨�Х����ȥǥ��쥯�ȥ꤬�������ĥ꡼�ǥ��쥯�ȥ�Ǥ�����, 
�ʲ��Υ��ޥ�ɤ��ǤĤȷ׻����¹Ԥ���ޤ�. 
�׻��ˤ� 10 ʬ���٤�����Ȼפ��ޤ�. 


  $ ./bin/arare -N=arare-mmc.conf


�׻���λ��, gpview �ʤɤ������Ԥʤ��ޤ�. 
�㤨�� gpview �����褹����, �ʲ��Τ褦�˥��ޥ�ɤ��Ǥ��ޤ�. 

  $ gpview --wsn 4 --range -4.0:4.0 MarsCond_PotTemp.nc@PotTempDist,t=1000

gpview �λ�����ˡ��
((<GPhys ���ޥ�ɥ��塼�ȥꥢ��|URL:http://ruby.gfd-dennou.org/tutorial/gpcommands/>))
�ʤɤ򻲾Ȥ��Ƥ�������.

MarsCond_VelZ.nc ���ѿ� VelZ �� 
MarsCond_DensCloud.nc ���ѿ� DensCloud 
MarsCond_PotTemp.nc ���ѿ� PotTempDist 
�˴ؤ��������Ԥʤ��Ȱʲ��Τ褦�ʿޤ������ޤ�. 


=end JA
=begin EN
== Overview

If you perform a calculation under default configuration, 
you have only to see this section "Overview". 
If not so, see the section "Details". 


If you finished compling the source for moist convection with major component 
condensation, you can execute the experiment. 
Note that command "make mmc" when compiling.
In order to execute the experiment, command as follows. 
It will take ten minutes or so. 

  $ ./bin/arare -N=arare-mmc.conf


After finishing the calculation, next step is visualization. 
For example, if you use gpview, enter a command as follows.

  $ gpview --wsn 4 --range -4.0:4.0 MarsCond_PotTemp.nc@PotTempDist,t=1000

If you want to know How to use gpview, for example, 
See
((<this page|URL:http://ruby.gfd-dennou.org/tutorial/gpcommands/>)). 

After drawing, you can get figures as follows. 

=end
=begin HTML
<tr>
<img src="mars-moistconv-denscloud.png" width="30%">
<img src="mars-moistconv-velz.png" width="30%">
<img src="mars-moistconv-pottemp.png" width="30%">
</tr>
=end HTML
=begin


=end EN

=begin JA

== �ܺ�

�ǥե���ȤȤϰۤʤ�����Ǽ¸���¹Ԥ�����, 
�ʲ��� 4 �ĤΥ��ƥåפǷ׻���Ԥ��ޤ�.

  * ���Ѥ���ץ��������򡦥��ԡ�
  * ����ѥ�᡼�����ѹ�
  * �������Υ���ѥ���
  * �׻��μ¹�

=== ���Ѥ���ץ��������򡦥��ԡ�

�����ޥ�徺�¸���Ԥʤ��٤ˤ�, 
���˼����褦�� src-mmc �ʲ��ˤ���ץ������ư�ǥ��ԡ����Ʋ�����. 

  $ cp main/arare_091012_masstest6-8_advlong-filter.f90 main/arare.f90
  $ cp io/historyfileio_mmconv2_masstest7.f90 io/historyfileio.f90
  $ cp physocs/densitycloud_turb-masstest2.f90 physics/densitycloud.f90
  $ cp physics/latentheat_all.f90 physics/latentheat.f90
  $ cp physics/masscondense_threshold-masstest4.f90 physics/masscondense.f90
  $ cp physics/adiation_balance4.f90 physics/radiation.f90
  $ cp setup/storepottemp_mmconv2.f90 setup/storepottemp.f90
  $ cp setup/storedenscloud_2d-5.f90 setup/storedenscloud.f90
  $ cp moist/moistbuoyancy_org.f90 moist/moistbuoyancy.f90 
  $ cp env/basicenv_smooth.f90 env/basicenv.f90
  $ cp numdiffusion_smooth_alpha1e-4.f90 util/numdiffusion.f90


=== ����ѥ�᡼�����ѹ�

�������ĥ꡼�ǥ��쥯�ȥ���β�������ʬ�ŷ���ή�Ѥ� namelist �ե����� 
arare-mmc.conf ���Խ����ޤ�. 

namelist �ե�����γƹ��ܤξܺ٤ˤĤ��Ƥ�
((<������|URL:../../doc/tutorial/namelist.htm>)) ������������. 

�Խ����ܤ�¿��¸�ߤ��ޤ���, 
�ʲ��Ǥ�(�Խ����Ƥ����꤬������ʤ��ȴ��Ԥ����)
���Τ����Ĥ���ԥå����åפ��ޤ�. 

* ��ʬ����

  TimeInt      = 1800.0d0

�ǥե���ȤǤ� 1800 sec �ȤʤäƤ��ޤ���, 
�׻�������ǽ�䵤ʬ�˹�碌��Ĺ��������û�������ꤷ�ƤߤƲ�����. 

* �ǡ������ϻ��ֳִ�

  TimeDisp     = 100.0d0

�ǥե���ȤǤ� 100 sec �ȤʤäƤ��ޤ���, 
�׻�������ǽ�䵤ʬ�˹�碌��Ĺ��������û�������ꤷ�ƤߤƲ�����. 

* �����ޥ�ο���

  DelMax       = 4.0d0

�ͤ򾮤������᤮���, ������������˥����ޥ뤬�Ĥ֤�Ƥ��ޤ��ޤ�. 
�ޤ��ͤ��礭�����᤮���, ʪ��Ū�˰�̣�Τ���������ʤ��ʤ�ޤ�. 
���ܾ���θ�����, 10 K ��Ķ���ʤ��褦�ˤ���Τ��ɤ��Ǥ��礦. 
�ޤ� DelMax �����ͤˤ�����䤿�����������ߤ���¸���Ԥʤ����Ȥ��Ǥ��ޤ�. 

* ����ˤ����륵���ޥ���濴����(��ľ����)���ΰ���Ф�����

  ZcRate       = 0.0d0

0 (��ü)���� 1 (��ü)�δ֤ǻ��ꤷ�ޤ�. 
�ǥե���ȤǤϥ�ǥ벼ü�˥����ޥ�����֤��Ƥ��ޤ�. 
�㤨�Ф��礦���濴�����֤������, 0.5d0 �Ȥ��ޤ�. 

* �׳�˰����

  SatRatioCr = 1.0d0

�ŷ뤬�Ϥޤ�ޤǤˤɤ�ۤɲ���Ѥˤʤ�ʤ���Фʤ�ʤ����򼨤��ѥ�᡼���Ǥ�. 
ʪ��Ū�˰�̣������׻���Ԥʤ��٤ˤ� 1 �ʾ�����ꤹ��ɬ�פ�����ޤ�. 
���ʤߤ˲����絤�Ǥ� 1.0 -- 1.35 �δ֤��ͤ�Ȥ�ȸ����Ƥ��ޤ�. 


=== �������Υ���ѥ���

�������ĥ꡼�ǥ��쥯�ȥ��, make mmc ���ޥ�ɤ�¹Ԥ��Ƥ�������. 
���˥���ѥ����¹Ԥ��Ƥ�����ˤ�, 
make mmc ��Ԥʤ����� make clean ��¹Ԥ��Ʋ�����. 

  $ make clean
  $ make mmc


=== �׻��μ¹�

�Ǹ�˷׻���¹Ԥ��ޤ�. 
�ꥹ�����ȥǡ���, ��ɽ�̥ꥹ�����ȥǡ����Ȥ����Ĥ��Υҥ��ȥ�ǡ���
�ե����뤬���Ϥ���ޤ�. 

�¹ԥץ���� bin/arare �ϰ����˲�����ꤷ�ʤ���, ��ȥǥ��쥯�ȥ����
arare.conf �� namelist �ե�����Ȥ����ɤ߹��ߤޤ�. 
��������ʬ�ŷ�׻���Ԥʤ����, arare-mmc.conf �� arare.conf �Ȥ��ƥ��ԡ�
���Ƥ���, �¹Ԥ��ޤ�. 

  $ cp arare-mmc.conf arare.conf
  $ ./bin/arare

�ޤ��� namelist �ե��������ꤷ�Ƽ¹Ԥ��뤳�Ȥ�Ǥ��ޤ�. 

  $ ./bin/arare -N=arare-mmc.conf

�¹Ի��ˤϰʲ��Υ��ץ�����Ϳ���뤳�Ȥ��Ǥ��ޤ�.

  -N=(NAMELIST �ե�����Υѥ�) �ޤ��� --namelist=(NAMELIST �ե�����Υѥ�)
     NAMELIST �ե�������ۤ˻��ꤹ��.

  -D �ޤ��� --debug
     �ǥХå���å���������Ϥ���.

  -H �ޤ��� --help
     �إ�ץ�å�������ɽ�����ƽ�λ����.

�ޤ�, ɸ����Ϥ����ɸ�२�顼���Ϥ���¸���뤳�Ȥ�Ǥ��ޤ�. 
�㤨�Х����뤬 bash �ξ��,

  $ /work/deepconv/bin/arare > arare.log 2> arare_error.log

�Ȥ��ޤ�. 
����ˤ��, ɸ����ϥե����� arare.log, ɸ�२�顼���ϥե����� 
arare_error.log ����������ޤ�. 


�׻�����λ�����, �ʲ��Τ褦�ʥǡ��������Ϥ���ޤ�. 


  MarsCond_restart2.nc     �ꥹ�����ȥե�����
  MarsCond_BasicZ.nc       ���ܾ��ѿ�
  MarsCond_DensCloud.nc    ��̩��
  MarsCond_Exner.nc        ̵�������ϴؿ�
  MarsCond_H2O-g.nc        �����������(�ܷ׻��Ǥϥ������Ϥ�������)
  MarsCond_H2O-s-Cloud.nc  ���庮����(�ܷ׻��Ǥϥ������Ϥ�������)
  MarsCond_H2O-s-Rain.nc   ���庮����(�ܷ׻��Ǥϥ������Ϥ�������)
  MarsCond_Kh.nc           ��ή�Ȼ�����(��ư��)
  MarsCond_Km.nc           ��ή�Ȼ�����(�����顼��)
  MarsCond_PotTemp.nc      ����
  MarsCond_SatRatio.nc     ˰����
  MarsCond_VelX.nc         ®��(X ��ʬ)
  MarsCond_VelZ.nc         ®��(Z ��ʬ)
  MarsCond_Zprof.nc        ������


����줿 netcdf �ե�����˴ؤ���, �����Ԥʤ���, 
�ֳ��ספǼ������褦�ʿޤ������ޤ�. 


=== �ꥹ�����ȥǡ�������μ¹�

���椫��³���Ʒ׻���Ԥʤ��������ˤ�, arare-mmc.conf ���Խ����ޤ�. 

�㤨�� MarsCond_restart2.nc ���˥ꥹ�����ȷ׻���Ԥʤ���������, 
�ʲ��Τ褦�˵��Ҥ��ޤ�. 

  (��)  InitFile    = ""   !����ͥե�����
  (��)  InitFile    = "MarsCond_restart2.nc"

�ꥹ�����ȷ׻���˽��Ϥ����ꥹ�����ȥե������̾�����Խ�����ˤ�, 
�ʲ��ι��ܤ��Խ����ޤ�. 
�ä�, �ꥹ�����Ȥ򳫻Ϥ���ե�����ȥꥹ�����ȷ׻���˽��Ϥ����ե�����
��̾����Ʊ��Ǥ�����, ��񤭤���ƥꥹ�����Ȥ򳫻Ϥ���ե����뤬̵����
�äƤ��ޤ��Τ���դ��Ʋ�����. 

  ReStartFile = "MarsCond_restart3.nc"

arare-mmc.conf ���Խ�������ä���, �ʲ��Υ��ޥ�ɤǼ¹Ԥ�ԤäƤ�������.
�ʤ�, ���Ϥ����ǡ����ե������̾���ϥǥե���ȤǤϷ�ޤäƤ���Τ�, 
�����׻������ǡ����ϥ�͡��ह�뤫, �̤Υǥ��쥯�ȥ�˰ܤ��Ƥ������Ȥ�
�����ᤷ�ޤ�. 

  $ ./bin/arare -N=arare-mmc.conf


=end JA

=begin EN

== Details

If you perform a calculation under non-default configuration, 
the experiment is performed with the following 4 steps:

  * Choosing and copying the programs
  * Changing the configuration
  * Compling
  * Execution of experiments

=== Choosing and copying the programs

For a Numerical experiment of ascending hot plume, 
you have to copy following programs in src-mmc directory.

  $ cp main/arare_091012_masstest6-8_advlong-filter.f90 main/arare.f90
  $ cp io/historyfileio_mmconv2_masstest7.f90 io/historyfileio.f90
  $ cp physocs/densitycloud_turb-masstest2.f90 physics/densitycloud.f90
  $ cp physics/latentheat_all.f90 physics/latentheat.f90
  $ cp physics/masscondense_threshold-masstest4.f90 physics/masscondense.f90
  $ cp physics/adiation_balance4.f90 physics/radiation.f90
  $ cp setup/storepottemp_mmconv2.f90 setup/storepottemp.f90
  $ cp setup/storedenscloud_2d-5.f90 setup/storedenscloud.f90
  $ cp moist/moistbuoyancy_org.f90 moist/moistbuoyancy.f90 
  $ cp env/basicenv_smooth.f90 env/basicenv.f90
  $ cp numdiffusion_smooth_alpha1e-4.f90 util/numdiffusion.f90


=== Changing the configuration

You may edit a namelist file for moist convection with major component 
condensation, "arare-mmc.conf". 

For details of namelist files, see
((<this page|URL:../../doc/tutorial/namelist.htm>)).

For example, you can edit the following entries easily. 

* Integration time

  TimeInt      = 1800.0d0

* Output time interval

  TimeDisp     = 100.0d0

* Amplitude of potential temperature of a plume

  DelMax       = 4.0d0

* Ratio of the altitude of a plume to vertical domain size

  ZcRate       = 0.0d0

* Critical saturation ratio

  SatRatioCr = 1.0d0


=== Compling

Before "make mmc", you should command "make clean".

  $ make clean
  $ make mmc



=== Execution of experiments

In order to execute the experiment, command as follows. 


If you give no argument, the execution program bin/arare read arare.conf
as a namelist file. 
When you perform experiments for Martian atmospheric convection with
major component condensation, you have to copy arare-mmc.conf as arare.conf, 
and execute bin/arare.

  $ cp arare-mmc.conf arare.conf
  $ ./bin/arare

Alternatively, you can also execute bon/arare by specifying a namelist file.

  $ ./bin/arare -N=arare-mmc.conf

When executing bin/arare, you can give following arguments.

  -N=(Path to namelist file) or --namelist=(Path to namelist file)
     Specifying a namelist file explicitly.

  -D or --debug
     Output the debug message.

  -H or --help
     Output the help message. 


Furthermore, you can save standard output and standard error output. 
If your shell is bash, command as follows. 

  $ ./bin/arare > arare.log 2> arare_error.log

Then you can get a standard output file "arare.log" and
a standard error output file "arare_error.log".


If you have finished calculation, you can obtain following data files. 

  MarsCond_restart2.nc     Restart file
  MarsCond_BasicZ.nc       Basic states
  MarsCond_DensCloud.nc    Cloud of density
  MarsCond_Exner.nc        The Exner function
  MarsCond_H2O-g.nc        Mixing ration for H2O vapor
  MarsCond_H2O-s-Cloud.nc  Mixing ration for H2O cloud
  MarsCond_H2O-s-Rain.nc   Mixing ration for H2O rain
  MarsCond_Kh.nc           Turbulent diffusion coefficient(momentum)
  MarsCond_Km.nc           Turbulent diffusion coefficient(scalar)
  MarsCond_PotTemp.nc      Potential Temperature
  MarsCond_SatRatio.nc     Saturation ratio
  MarsCond_VelX.nc         Horizontal velocity
  MarsCond_VelZ.nc         Vertical velocity
  MarsCond_Zprof.nc        Diagnosis

You can get figures shown above by vizualisation for data files. 

=end EN


=begin JA
== ����ʸ��
=end JA

=begin EN
== References
=end EN

=begin

* M. Odaka, T. Kitamori, K. Sugiyama, K. Nakajima and Y.-Y. Hayashi, 2006:
  Numerical simulation of Martian atmospheric convection, 
  Proc. of the 20th ISAS Atmospheric Science Symposium, JAXA/ISAS, 
  103 -- 106

* Klemp, J. B., Wilhelmson, R. B., 1978: 
  The simulation of three-dimensional convective storm dynamics", 
  J. Atmos. Sci., 35, pp. 1070 -- 1096.

* Tobie, G., Forget, F., Lott, F., 2003: 
  Numerical simulation of winter polar wave clouds observed by
  Mars Global Surveyor Mars Orbiter Laser Altimeter", 
  Icarus, 35, 33 -- 49.

=end

=begin HTML
<hr />
<small>
  $Id: mars-plumetest.rd,v 1.1 2011/06/21 15:42:02 sugiyama Exp $
</small>
=end HTML

