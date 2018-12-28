=begin JA

= ��ñ�ʲ��ϡ��Ļ벽

# * ���� ���� (morikawa), Ǽ¿ ů�� (noda), �ⶶ ˧�� (yot), �ݹ� ���� (takepiro)
# * ���� �̰�ϯ (sugiyama) [copy from dcpam]
#   * $Id: visualization.rd,v 1.2 2014/03/04 14:02:28 sugiyama Exp $

=end JA
=begin EN

= First step analysis and visualization

# * Yasuhiro MORIKAWA (morikawa), Satoshi NODA (noda), Yoshiyuki * O. Takahashi (yot) Shin-ichi Takehiro (takepiro)
# * Ko-ichiro SUGIYAMA (sugiyama) [copy from dcpam]
#   * $Id: visualization.rd,v 1.2 2014/03/04 14:02:28 sugiyama Exp $

=end EN


=begin JA
== ���ϡ��Ļ벽�ġ���ν���

deepconv/arare5 �Ǥ������Ϥ���ե�����Ȥ���
((<Gtool4 NetCDF ����|URL:http://www.gfd-dennou.org/library/gtool>))
�˴�Ť��� NetCDF �ǡ����򰷤��ޤ�. 

���ͼ¸��η�̤���ϡ��Ļ벽���뤿��ˤ�, NetCDF
�ǡ������갷�����ȤΤǤ�����ϡ��Ļ벽�ġ��뤬ɬ�פǤ�. 
�����Ǥ�, ((<��Ǿ Ruby �ץ�������|URL:http://ruby.gfd-dennou.org/index-j.htm>))
�����󶡤���� ((<Gphys|URL:http://www.gfd-dennou.org/library/ruby/products/gphys/>))
��Ȥä��Ļ벽�����Ҳ𤷤ޤ�. 

=== �Ļ벽�ġ���Υ��󥹥ȡ���

((<��ǾRuby������ ���󥹥ȡ��륬����|URL:http://www.gfd-dennou.org/arch/ruby/tutorial/install/index-j.html>))
�򻲾Ȥ��Ƥ�������. 

=end JA

=begin EN
== Preparation of analysis and visualization tools

"deepconv/arare5" input/output NetCDF data based on 
((<Gtool4 NetCDF Conventions|URL:http://www.gfd-dennou.org/library/gtool/index.htm.en>))

Analysis and visualization tools for NetCDF data are
need in order to analyze and visualize results of numerical experiments.
Here, ((<Gphys|URL:http://www.gfd-dennou.org/library/ruby/products/gphys/>))
provided from 
((<Dennou Ruby Project|URL:http://ruby.gfd-dennou.org/index.htm>))
is used via irb. 
Details of usage can be found ((<here|URL:http://www.gfd-dennou.org/library/ruby/tutorial/index-e.html>)).

=== Installation of tool

See ((<Dennou Ruby Products Installation Guide|URL:http://www.gfd-dennou.org/arch/ruby/tutorial/install/>)). 

=end EN

=begin JA

== GPhys/GGraph �ˤ����ϤȲĻ벽

�����Ǥ�, 
((<����ή�μ¸�|URL:./exp-s93.htm>))
������줿�ǡ����� GPhys/GGraph ���Ѥ��ƲĻ벽���Ƥߤ뤳�Ȥˤ��ޤ�. 

�ޤ� irb ��ư���Ƥ�������. 

  $ irb

�ʲ��Τ褦�� irb �Υץ��ץȤ�ɽ������ޤ�. 

  irb(main):001:0>

���Υץ��ץȤ�, �ʲ��Τ褦�˥��ޥ�ɤ��Ǥ��ޤ�. 
��ü�ο����Ϲ��ֹ��, �Ǥ�ɬ�פϤ���ޤ���.

  1: require "numru/ggraph"
  2: include NumRu
  3: gphys = GPhys::IO.open('denscurrent-dry_PTemp.nc', 'PTemp').cut('y'=>0.0)
  4: DCL.gropn(4)
  5: DCL.sgpset('lcntl', false) ; DCL.uzfact(0.7)
  6: GGraph.tone gphys

irb �Υץ��ץȤˤ����� quit ���ǤĤ� irb ��λ���뤳�Ȥ��Ǥ��ޤ�. 

�����Ǥ� Temp.nc �Ȥ����ե��������� Temp �Ȥ����ѿ����ɤ߹���, 
�޼���ԤäƤ��ޤ�. 

((<"IMG:exp-s93_img02.png">))

PTemp �� x, y, z, t (����) �� 4 �����ǡ����Ǥ���, 
���μ¸��Ͽ�ʿ��ľ 2 �����ǹԤäƤ���Τ�, 
3 ���ܤ� y = 0 �Ȥ��Ƥ��ޤ�. 
������ꤷ�ʤ��ȺǸ�μ����˴ؤ��Ƥϼ�ưŪ�� 1 ���ܤ����Ǥ����򤵤�ޤ�. 
�������äƤ��οޤ� t=0 �Ǥο�ʿ��ľ 2 �����Ǥβ��̤򼨤��Ƥ��뤳�Ȥˤʤ�ޤ�

�ޤ�, ³����, ���Τ褦�˻������ꤹ�뤳�Ȥ�, 
�ۤʤ����Ǥκǲ��ؤβ���ʬ�ۤ��������Ȥ��Ǥ��ޤ�. 

  7: GGraph.tone gphys.cut('t'=>900)

((<"IMG:exp-s93_img03.png">))

��λ������ˤ� 

  8: DCL.grcls
  9: quit

�Ȥ��ޤ��礦. 

�������������Ǥʤ�, ���Ϥ�Ԥ����Ȥ�Ǥ��ޤ�. 
��Ȥ��Ʋ��̤ȥ������ʡ��ؿ����鲹�٤�׻���, �޼����Ƥߤޤ��礦. 

  1: require "numru/ggraph"
  2: include NumRu
  3: gphys1 = GPhys::IO.open('denscurrent-dry_PTempAll.nc', 'PTempAll').cut('y'=>0.0)
  4: gphys2 = GPhys::IO.open('denscurrent-dry_ExnerAll.nc', 'ExnerAll').cut('y'=>0.0)
  5: temp = gphys1 * gphys2
  6: DCL.gropn(1)
  7: DCL.sgpset('lcntl', false) ; DCL.uzfact(0.7)
  8: GGraph.tone temp.cut('t'=>900,'x'=>25e3..45e3)
  9: GGraph.contour temp.cut('t'=>900,'x'=>25e3..45e3),false
  10: GGraph.color_bar

((<"IMG:exp-s93_img04.png">))

�� 5 ���ܤǲ��٤η׻��򲹰̤ȥ������ʡ��ؿ�����ԤäƤ��ޤ�. 
�޼������Ƥ���Τ� x = 25 ~ 45 km �Ǥ�. 

���Τ褦��, GPhys/GGraph ���Ѥ����¿�̤ʲ��ϤȲĻ벽��¸��Ǥ��ޤ�. 
������Ĺ���ʤäƤ����� irb �ǥ��󥿥饯�ƥ��֤˹Ԥ�������, 
���ǥ������Ѥ��ƥ�����ץȥե�������������ΨŪ�Ǥ�������Ѥ��ưפˤʤ�ޤ�. 
�����٤ʲ��ϡ��Ļ벽��Ԥ��ݤˤ�, 
((<GPhys ���塼�ȥꥢ��|URL:http://ruby.gfd-dennou.org/products/gphys/tutorial/>))
�򻲾Ȥ��Ƥ�������. 

=end JA

=begin EN

== Analysis and visualization with GPhys/GGraph

Under construction

=end EN

=begin JA

== GPhys/gp���ޥ�ɤˤ��Ļ벽

�����Ǥ�, 
((<Polvani et al. (2004) �η����԰�����ư�¸�|URL:./exp-p04.htm>))
������줿�ǡ����� GPhys ��°�� gp ���ޥ�ɤ��Ѥ��ƲĻ벽���Ƥߤ뤳�Ȥˤ��ޤ�. 
���٤Υǡ������ɤ߼��޼�����ˤ�, 

  % gpview denscurrent-dry_PTemp.nc@PTemp,y=0

�����Ϥ��ޤ�. ����� denscurrent-dry_PTemp.nc �Ȥ����ե��������� 
PTemp �Ȥ����ѿ����ɤ߹���, �޼�����Ȥ������ޥ�ɤǤ�. 

((<"IMG:exp-s93_img05.png">))

PTemp �� x, y, z, t (����) �� 4 �����ǡ����Ǥ���, 
���μ¸��Ͽ�ʿ��ľ 2 �����ǹԤäƤ���Τ�, y = 0 �Ȥ��Ƥ��ޤ�. 
������ꤷ�ʤ��ȺǸ�μ����˴ؤ��Ƥϼ�ưŪ�� 1 ���ܤ����Ǥ����򤵤�ޤ�. 
�������äƤ��οޤ� t=0 �Ǥο�ʿ��ľ 2 �����Ǥβ��̤򼨤��Ƥ��뤳�Ȥˤʤ�ޤ�

���̤��Ѥ���, ���� t=900, x = 25 ~ 50 km �ˤ��������, 

  % gpview denscurrent-dry_PTemp.nc@PTemp,y=0,x=25e3:50e3,t=900 --nocont 

��, ����ޤǶ��ڤä����̤���ꤹ�뤳�Ȥ��Ǥ��ޤ�. 

((<"IMG:exp-s93_img06.png">))

���˥᡼�������ñ�˸��뤳�Ȥ��Ǥ��ޤ�. 

  % gpview denscurrent-dry_PTemp.nc@PTemp,y=0,x=25e3:50e3 --anim t --nocont

��ɸ�� 3 �Ļ��ꤷ�� 1 �����ǡ����ˤ�����ޤ�������դ������ޤ�. 

  % gpview denscurrent-dry_PTemp.nc@PTemp,y=0,x=25e3:50e3,t=900 --exch

--exch ���ץ����Ͻļ��Ȳ��������촹��������ؼ������ΤǤ�. 

((<"IMG:exp-s93_img07.png">))

gpview �ˤ�¾�ˤ⤤����ʥ��ץ���󤬤���ޤ�. 
gpview --help �Ȥ���ȥ��ץ����ȻȤ������㤬ɽ������ޤ�. 
gp ���ޥ�ɥ��꡼���ˤ�¾�ˤ⤤����ʤ�Τ��Ѱդ���Ƥ��ޤ�. 
��ʤ�Τϰʲ����̤�Ǥ�. 

: gpvect
  2 �����٥��ȥ�ޤ�ɽ��

: gpprint
  �ǡ����ο��ͤν���ɽ��

: gplist
  �ե�����˳�Ǽ����Ƥ����ѿ��Υꥹ�Ȥ�ɽ��

: gpmaxmin
  �ǡ����κ��硦�Ǿ��ͤ�ɽ��

gp ���ޥ�ɥ��꡼���� 1 �����ϤǤ����˷�̤�ɽ���Ǥ���Τ���ħ�Ǥ�. 
���Τ��᥯���å���å���׻�������ǤΥǡ��������å����������Ǥ�. 
�������ʤ���ʣ�����ѿ����Ȥ߹�碌�����Ϥ�Ļ벽�ϤǤ��ޤ���. 
�ܳ�Ū�ʥǡ������ϤȲĻ벽��Ԥ��ˤ�, 
((<"GPhys/GGraph �ˤ����ϤȲĻ벽">))��Ŭ���Ƥ���Ǥ��礦. 

=end JA

=begin EN

== Visualization with GPhys/gpcommands

Under construction

=end EN

=begin HTML
<hr />
<small>
  $Id: visualization.rd,v 1.2 2014/03/04 14:02:28 sugiyama Exp $
</small>
=end HTML

