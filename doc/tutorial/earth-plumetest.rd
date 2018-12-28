=begin TOPLINK
[((<English|URL:mars-plumetest.htm.en>)) |
((<Japanese|URL:mars-plumetest.htm>))]
[((<GFD Dennou Club|URL:http://www.gfd-dennou.org>)) |
((<deepconv Project|URL:http://www.gfd-dennou.org/library/deepconv>))]
=end TOPLINK


=begin JA

= �ϵ��絤��ﲼ�ǤΥ����ޥ�徺�¸�

#* ����
#  * 2011/03/04(�����̰�ϯ) ��������


# * ���� ã�� (yamasita), ���� �̰�ϯ (sugiyama), ���� ���� (odakker)
#   * $Id: earth-plumetest.rd,v 1.1 2011/06/21 15:42:01 sugiyama Exp $

=end JA
=begin EN

= A Numerical experiment of ascending hot plume

#* History
#  * 2011/03/01 (Tatsuya Yamashita) Initial release


# * Tatsuya YAMASHITA (yamasita), Ko-ichiro SUGIYAMA (sugiyama), Masatsugu ODAKA (odakker)
#   * $Id: earth-plumetest.rd,v 1.1 2011/06/21 15:42:01 sugiyama Exp $

=end EN


=begin JA

�����Ǥ�, ñ�ȤΥۥåȥץ�塼���徺������¸���¹Ԥ�����ˡ���������ޤ�. 

���η׻��ˤϰʲ���ʪ���������Ѥ��Ƥ��ޤ�. 

  * ��ή�Ȼ�
    * 1.5 ���Υ������㡼��ǥ� (Klemp and Wilhelmson, 1978)
  * �ŷ롦��ȯ
    * Kessler �Υѥ�᥿�ꥼ�������

=end JA
=begin EN

A method to perform an test experiment
of an ascending hot plume is described. 

Following physical processes are used in this experiment.

  * Turbulent mixing
    * 1.5 order closure (Klemp and Wilhelmson, 1978)
  * Condensation and evaporation
    * Kessler parametarization

=end EN


=begin JA
== ����

�ܼ¸��ϥǥե���Ȥ�����ǹԤʤ����Ȥ�����ޤ�. 
�ǥե��������Ǥη׻��򤹤���Ϥ��Ρֳ��ספΤߤ��ɤ�н�ʬ�Ǥ�. 
�ǥե���ȤȤϰۤʤ�����Ƿ׻��򤹤����, �־ܺ١פ�����������. 

�㤨�Х����ȥǥ��쥯�ȥ꤬�������ĥ꡼�ǥ��쥯�ȥ�Ǥ�����, 
�ʲ��Υ��ޥ�ɤ��ǤĤȷ׻����¹Ԥ���ޤ�. 
�׻��ˤ� 10 ʬ���٤�����Ȼפ��ޤ�. 


  $ ./bin/arare -N=arare-earth.conf


�׻���λ��, gpview �ʤɤ������Ԥʤ��ޤ�. 
�㤨�� gpview �����褹����, �ʲ��Τ褦�˥��ޥ�ɤ��Ǥ��ޤ�. 

  $ gpview --wsn 4 arare-earth_PotTemp.nc@PotTemp --anim t

gpview �λ�����ˡ��
((<GPhys ���ޥ�ɥ��塼�ȥꥢ��|URL:http://ruby.gfd-dennou.org/tutorial/gpcommands/>))
�ʤɤ򻲾Ȥ��Ƥ�������.


=end JA
=begin EN
== Overview

If you perform a calculation under default configuration, 
you have only to see this section "Overview". 
If not so, see the section "Details". 


If you finished compling the source for moist convection with major component 
condensation, you can execute the experiment. 
In order to execute the experiment, command as follows. 
It will take ten minutes or so. 

  $ ./bin/arare -N=arare-earth.conf


After finishing the calculation, next step is visualization. 
For example, if you use gpview, enter a command as follows.

  $ gpview --wsn 4 arare-earth_PotTemp.nc@PotTemp --anim t

If you want to know How to use gpview, for example, 
See
((<this page|URL:http://ruby.gfd-dennou.org/tutorial/gpcommands/>)). 

=end EN

#=begin HTML
#<tr>
#<img src="mars-moistconv-denscloud.png" width="30%">
#<img src="mars-moistconv-velz.png" width="30%">
#<img src="mars-moistconv-pottemp.png" width="30%">
#</tr>
#=end HTML


=begin JA
== ����ʸ��
=end JA

=begin EN
== References
=end EN

=begin

* Klemp, J. B., Wilhelmson, R. B., 1978: 
  The simulation of three-dimensional convective storm dynamics", 
  J. Atmos. Sci., 35, pp. 1070 -- 1096.


=end

=begin HTML
<hr />
<small>
  $Id: earth-plumetest.rd,v 1.1 2011/06/21 15:42:01 sugiyama Exp $
</small>
=end HTML

