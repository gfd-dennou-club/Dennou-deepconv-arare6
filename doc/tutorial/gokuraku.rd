#= Totorial of deepconv/arare: rakuraku deepconv
#
# Authors::   SUGIYAMA Ko-ichiro, ODAKA Masatsugu
# Version::   $Id: gokuraku.rd,v 1.6 2014/02/26 07:13:10 sugiyama Exp $
# Tag Name::  $Name:  $
# Copyright:: Copyright (C) GFD Dennou Club, 2006. All rights reserved.
# License::   See COPYRIGHT[link:../../COPYRIGHT]
#
#
########################################################################
#
=begin TOPLINK
[((<GFD Dennou Club|URL:http://www.gfd-dennou.org>)) |
((<deepconv Project|URL:http://www.gfd-dennou.org/library/deepconv>))]
[((<deepconv Reference Manual|URL:../index.htm>))]
=end TOPLINK

=begin

= �����餯 deepconv/arare5
#= deepconv/arare5 ���Ѥμ��
#* ���� �̰�ϯ, ���� ����, ���� ã��
#  * 2014/02/05  (���� �̰�ϯ) ����
#  * 2012/05/09  (���� ����) ����
#  * 2011/03/03  (���� ã��) ����
#  * 2011/02/24  (���� ã��) ����
#  * 2009/03/06  (���� ã��) ����
#  * 2008/06/18  (���� ����) ����
#  * 2006/11/20  (���� ����) ����
#  * 2006/10/18  (���� ����) ����
#  * 2006/09/12  (���� ����) ��������

== �Ϥ����

����ʸ��� deeoconv/arare ���Ѥ��Ƽ�ڤ˼¸���Ԥ�����Υ��塼�ȥꥢ��
�Ǥ�. 

#�¹ԥץ����ϥ������ĥ꡼�ǥ��쥯�ȥ� (�����Ǥ� /work/deepconv 
#��Ÿ������Ƥ���Ȥ��ޤ�) ľ���� bin �ǥ��쥯�ȥ�ʲ��ˤ���Ȥ��ޤ�. 

== deepconv/arare5 �Υӥ��

((<deepconv ���󥹥ȡ���μ��|URL:../../INSTALL.htm>)) �򻲹ͤ�, 
deepconv/arare5 �Υӥ�ɤ�ԤäƤ�������. 
�֥��󥹥ȡ���μ��פΡ֥������Υ���ѥ���פޤǹԤäƤ�������. 
((<�����Ĥ��Υ���ѥ���˴ؤ�����ս�|URL:./compiler_note.htm>)) �⻲�Ȳ�����. 

�ӥ�ɤ���λ�����, "src/main" �ǥ��쥯�ȥ�ʲ���, arare, arare_initdata �Ȥ��ä��¹ԥե����뤬��������ޤ�. 
�ޤ�, �����Ĥ��Υ���ץ� NAMELIST �ե����� (��ĥ�Ҥ� .conf �Υե�����) �� "exp_setup_files" ���Ѱդ���Ƥ��ޤ�.


== �¸��μ¹�

����, �����Ĥ��μ¸��μ¹���ˡ�ˤĤ��Ƥ������򵭤��ޤ�. �Ϥ�� deepconv ��Ȥ��ͤ�, �ޤ��ɤ줫 1 �Ĥμ¸� (�㤨�С�XXXXXXXXXXXXXX��) ��¹Բ�����. 

���ʤߤ�, �ɤμ¸���¹Ԥ���ˤϰʲ��� 4 �ĤΥ��ƥåפǹԤ��ޤ�.

* �¸��ǥ��쥯�ȥ�ν���
* ����ͤν���
* �¸��ѥǡ����ν���
* �¸��μ¹�

�ʤ�, �¸��ǥ��쥯�ȥ�ν�����, ɬ������ɬ�פ���ޤ���. �����Ǥ�, ����¸���������̤�¾�Τ�ΤȺ����äƤ��ޤ����Ȥ��ɤ������, �Ƽ¸����Ȥ˥ǥ��쥯�ȥ��������Ƥ��ޤ�.

===����ή�μ¸� (Straka et al., 1993)

����ή�μ¸���¹Ԥ�����ˡ��((<������|URL:./exp-s93.htm>))���������ޤ�. 

===����ӥ�إ��ۥ���԰���μ¸� 

����ӥ�إ��ۥ���԰���μ¸���¹Ԥ�����ˡ��((<������|URL:./catseye.htm>))���������ޤ�. 


#* �����ޥ�ξ徺 (�ŷ�ʤ�, ��ή: Klemp and Wilhelmson)
#
#   arare-thermal-dry_init-data.conf
#   arare-thermal-dry.conf

#* �����ޥ�ξ徺 (�ŷ�: Kessler, H2O �Τ�, ��ή: Klemp and Wilhelmson)
# 
#   arare-thermal-moist_init-data.conf
#   arare-thermal-moist.conf

#* �����ޥ�ξ徺 [������] (�ŷ�: Kessler, H2O, NH3, NH4SH, ��ή: Klemp and Wilhelmson)
#
#   arare-jupiter_init-data.conf
#   arare-jupiter.conf

#* �����ޥ�ξ徺 [������] (�ŷ�: �Ȼ���Ĺ, ��ή: Klemp and Wilhelmson)
#
#   arare-mars_init-data.conf
#   arare-mars.conf

#* Takemi (2007) �ν����
#
#   arare-takemi_init-data.conf

#* Yamasaki (1983) �ν����
#
#   arare-yamasaki_init-data.conf

#* ���С������ (deepconv/arare4) �Υҥ��ȥ꡼�ե����뤫����������
#
#   arare-from-arare4_init-data.conf 

#�ʲ�, �����Ĥ��μ¸��μ¹���ˡ�ˤĤ����������ޤ�. 
#
#=== ���������ޥ�徺�¸�
#
#�����絤���Ϥ�����ﲼ��, �����ޥ�徺�¸���¹Ԥ�����ˡ��
#((<������|URL:./jupiter-plumetest.htm>)) ���������ޤ�. 
#
#=== �ϵ奵���ޥ�徺�¸�
#
#�ϵ��絤���Ϥ�����ﲼ��, �����ޥ�徺�¸���¹Ԥ�����ˡ��
#((<������|URL:./earth-plumetest.htm>)) ���������ޤ�. 
#
#=== ���������ޥ�徺�¸�
#
#����¾ (2006) �ǹԤ�줿, ���������ޥ�徺�¸���¹Ԥ�����ˡ��
#((<������|URL:../../doc-mmc/tutorial/mars-plumetest.htm>)) ���������ޤ�. 


== ��ñ�ʲ��ϡ��Ļ벽

��ñ�ʲ��ϡ��Ļ벽�ˤĤ��Ƥ�, ((<��ñ���ϡ��Ļ벽|URL:./visualization.htm>)) �򻲹ͤˤ��Ƥ�������. 

#== �¸������ѹ�
##
#�¸������ѹ�����ˤ� namelist ���Խ�����ɬ�פ�����ޤ�. 
#namelist �ե�����γƹ��ܤˤĤ��Ƥ� 
#((<������|URL:namelist.htm>)) ���������ޤ�. 

== �¹ԥץ����ι���

src �ǥ��쥯�ȥ�ʲ����Խ���Ԥä����, 
�ʲ��Τ褦�ˤ��뤳�ȤǼ¸��ѥǥ��쥯�ȥ�����Ƥ򹹿����뤳�Ȥ���ǽ�Ǥ�. 

�ޤ�, �������ĥ꡼ľ����

 $ make clean

��¹Ԥ��ޤ�. 

����
((<���󥹥ȡ���μ��|URL:../../INSTALL.htm>)) ��
�֥��󥹥ȡ���μ��פΡ֥������Υ���ѥ���פ򻲹ͤ�, 
�Ԥʤ������׻��˱����ƥ������Υ���ѥ������ٹԤʤ��ޤ�. 

#== ����ʸ��
#
#* M. Odaka, T. Kitamori, K. Sugiyama, K. Nakajima and Y.-Y. Hayashi, 2006:
#  Numerical simulation of Martian atmospheric convection,
#  Proc. of the 20th ISAS Atmospheric Science Symposium, JAXA/ISAS,
#  103 -- 106


=end
