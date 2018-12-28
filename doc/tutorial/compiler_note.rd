=begin JA

= �����Ĥ��Υ���ѥ���˴ؤ�����ս�

# * ���� ����, �ⶶ ˧��
# * ���� �̰�ϯ, dcpam ���饳�ԡ�
#   * $Id: compiler_note.rd,v 1.2 2014/07/08 07:47:12 sugiyama Exp $

=end JA
=begin EN

= Notes about some compilers

# * Yasuhiro Morikawa, Yoshiyuki O. Takahashi
# * Ko-ichiro Sugiyama, copy from dcpam
#   * $Id: compiler_note.rd,v 1.2 2014/07/08 07:47:12 sugiyama Exp $

=end EN

=begin HTML
<a name="ifort">
<a name="g95">
=end HTML

=begin JA


== PGI Fortran

=== ����ѥ��륪�ץ����

��ư�������黻�����٤��ݾڤ��뤿��Υ��ץ�����ɬ���դ��Ʋ�����. 
deepconv �ˤ����Ƥ�, ����ѥ��顼�κ�Ŭ���򶯤�����������, 
�׻���̤����������ʤ��㤬��𤵤�Ƥ��ޤ�. 

   -Kieee


== Intel Fortran

=== ����ѥ��륪�ץ����

��ư�������黻�����٤��ݾڤ��뤿��Υ��ץ�����ɬ���դ��Ʋ�����. 
deepconv �ˤ����Ƥ�, ����ѥ��顼�κ�Ŭ���򶯤�����������, 
�׻���̤����������ʤ��㤬��𤵤�Ƥ��ޤ�. 

   -fp-model strict -prec-div
   (�⤷����, -fp-model precise -prec-div)



== Intel Fortran, G95 Fortran

=== �������ơ�����󥨥顼

OS ����ӥ����ͥ�ΥС������ˤ�äƤ�, �����٤�夲���׻���
�Ԥä��ݤ˥������ơ�����󥨥顼���������礬����ޤ�.

���ξ��ˤ�, ����Υ����å��ΰ�λ������̤κ����ͤ����Ū������
���ꤵ��Ƥ����ǽ��������ޤ�. 

¿���ξ��ˤϤ��λ������̤κ����ͤ�
�桼���Ǥ������ǽ�ʥ��եȥ�ߥåȤȸƤФ���ΤǤ��뤿��, 
�ʲ��Τ褦�˺����ͤ����䤹�褦�����ǽ�Ǥ�. 
������, �����ƥ�¦�Ƿ����Ƥ�������ͤǤ���ϡ��ɥ�ߥåȰʾ��
�礭�����եȥ�ߥåȤ����ꤹ�뤳�ȤϤǤ��ޤ���. 


==== �����å��ΰ�Υ��եȥ�ߥåȤ�ϡ��ɥ�ߥåȤޤ����䤹

: csh, tcsh �ǥ����å��ΰ�Υ��եȥ�ߥåȤ������ͤ��ѹ�����

  �����å��ΰ�Υϡ��ɥ�ߥåȤ�Ĵ�٤ޤ�
  
    > limit -h | grep stack
    stacksize    XXXXXXXX

  ������ XXXXXXXX �Ȥ���ɽ�������Τ�
  �����å��ΰ�Υϡ��ɥ�ߥåȤǤ�.
  (���ͤޤ��� "unlimited" ��ɽ������ޤ�). 
  ���� XXXXXXXX �򥹥��å��ΰ�Υ��եȥ�ߥåȤ����ꤷ�ޤ�.

    > limit stacksize XXXXXXXX

: sh, bash �ǥ����å��ΰ�Υ��եȥ�ߥåȤ������ͤ��ѹ�����

  �����å��ΰ�Υϡ��ɥ�ߥåȤ�Ĵ�٤ޤ�
  
    $ ulimit -Ha |grep stack
    stack size            (kbytes, -s) XXXXXXXX

  ������ XXXXXXXX �Ȥ���ɽ�������Τ�
  �����å��ΰ�Υϡ��ɥ�ߥåȤǤ�.
  (���ͤޤ��� "unlimited" ��ɽ������ޤ�). 
  ���� XXXXXXXX �򥹥��å��ΰ�Υ��եȥ�ߥåȤ����ꤷ�ޤ�.

    $ ulimit -s XXXXXXXX

##
## ��ƣ���󤫤� Mac �ξ��ˤĤ��ƾ�����餦����.
##

# �ʲ��ϥ��� (����������νŻ᰸�˽񤤤��᡼����ȴ��)
#
#== �ץ���बư���ʤ��ä�����
#
#�ץ���बư���ʤ��ä������ϥ����å��ΰ�˳�����Ƥ���������­
#�Ǥ�. �¥����������Ѥ�Ǥ��Ƥ�, �ƥ桼����˥���λ������̤˸³�
#�����ꤵ��Ƥ���, limit �� ulimit -a ���ޥ�ɤǤ��θ³���ɽ���Ǥ��ޤ�. 
#�ɤ���� Linux Kernel 2.6 �ϤǤϥ����å��˳�����Ƥ��륽�եȥ�ߥå�
#���� 8 MB �ˤʤäƤ��뤿��, ������­�ǥ������ơ�����󥨥顼��ȯ��
#���Ƥ����褦�Ǥ�.
#
#���ʤߤ�, ����ˤϥ����å��ΰ�ʳ��˥ҡ����ΰ�⤢��, �ץ���ब��
#���å��ΰ����Ѥ��뤫�ɤ����ϥץ����˽��� (�ä˳��ѿ������ʸ��
#�ɤΤ褦�˽񤯤�) �˵��ޤ�. dcpam5 �Ǥ�, �ƥ���ѥ��뤻���˲����٤�
#�Ѥ��뤿��, �ޤ���Ū�ʲ�������ͥ�褹��ʤɤ���ͳ���������ưŪ���դ˼�
#ư������Ѥ��Ƥ��ޤ�. ��������˴ؤ��Ƥϥ���Υ����å��ΰ褬���Ѥ���
#�뤳�Ȥˤʤ�ޤ�. �㤨�� T42L20 �ξ��, ���ѥ��ꥵ���������Τ��� 72
#MB �Ǥ�. ���Τ��������å��ΰ�����Ѥ���ʬ�����եȥ�ߥåȤ� 8 MB ��Ķ
#���뤿��˥��顼�������Ƥ��ޤ�.
#
#�ʤ�, ��ư����ʤɤ˴ؤ��Ƥ�, ����ģ Fortran ɸ�ॳ���ǥ��󥰥롼��
#http://www.mri-jma.go.jp/Project/mrinpd/coderule.html �Ρ�3.2 ưŪ���
#�դ��פ����ͤˤʤ뤫�Ȼפ��ޤ�. (dcpam �Ǥ����Ū�ˤϤ��Υ����ǥ���
#�롼����򤷤Ƥ��ޤ�).
#
#


�ޤ�, Intel Fortran �� ver.10.0 �ʾ�Ǥϥ���ѥ��륪�ץ���� (FFLAGS) ��

  -heap-arrays

���ɲä��Ƥ�������. ����ϼ��Τ褦�ʻ���ˤ��ޤ�.

Intel Fortran �� ver.10.0 �ʾ�Ǥϰ��Ū�ʥ��������Ƥ�
�����å��ؤγ���դ��뤳�Ȥ������ޤ���.
���Τ���, �嵭�Τ褦�� limit �ޤ��� ulimit �ˤ�����¤�ˤ�Ƥ�
�������ơ�����󥨥顼���Ф뤳�Ȥ�����ޤ�.
-heap-arrays ����ꤹ��Ȱ��Ū�ʥ��������Ƥ�
�ҡ����ΰ��ˤ����ʤ���褦�ˤʤ�Τ�, �������ơ�����󥨥顼
������Ǥ��ޤ�. 

���͡�http://software.intel.com/en-us/articles/intel-fortran-compiler-increased-stack-usage-of-80-or-higher-compilers-causes-segmentation-fault/


=end JA

=begin EN
== Intel Fortran, G95 Fortran

=== Segmentation fault

When a high resolution calculation is performed, a segmentation fault
might be caused according to version of OS and kernel.

In this case, there is a possibility that the maximum value of the use
capacity of a stack area of the memory is set comparatively small. 

In many cases, the maximum value of this use capacity is the one that
is called "soft limit" that users can set. Therefore, it can be set
that the maximum value is increased as follows.
However, a soft limit cannot be greatly set more than a "hard limit"
that is the maximum value that has been decided on the system side.

==== A soft limit of a stack area is increased to a hard limit

: With csh or tcsh, A setting value of a soft limit of a stack area is changed

  A hard limit of a stack area is examined as follows.
  
    > limit -h | grep stack
    stacksize    XXXXXXXX

  "XXXXXXXX" is a value of a hard limit.
  (Numerical value or "unlimited" is displayed).
  This "XXXXXXXX" is set to a soft limit of a stack area as follows. 

    > limit stacksize XXXXXXXX

: With sh or bash, A setting value of a soft limit of a stack area is changed

  A hard limit of a stack area is examined as follows.
  
    $ ulimit -Ha |grep stack
    stack size            (kbytes, -s) XXXXXXXX

  "XXXXXXXX" is a value of a hard limit.
  (Numerical value or "unlimited" is displayed).
  This "XXXXXXXX" is set to a soft limit of a stack area as follows. 

    $ ulimit -s XXXXXXXX


=end EN

=begin HTML
<a name="frt">
=end HTML

=begin JA
== Fujitsu Fortran

=== NAMELIST �ɤ߹��ߥ��顼

�ե�����˵��Ҥ���� NAMELIST ��¿������, �ޤ�˰�����
NAMELIST ���ɤ߹��ޤ�ʤ���礬����褦�Ǥ�. 

�ե�����ˤϤ����Ƚ񤭹��ޤ�Ƥ���Ϥ��ʤΤ�, NAMELIST �ξ��󤬼¹�
�ץ�����ȿ�Ǥ���ʤ����ˤ�, �¹ԥץ����ˤ�äƽ��Ϥ������
������������å����Ƥ�������. �⤷��ʲ��Τ褦�ʥ�å�������ɽ�������
�����, �����ǽҤ٤�褦�ʾɾ�������Ƥ��뤳�Ȥˤʤ�ޤ�.

 !*** WARNING [XxxxNmlRead] ***  NAMELIST group "xxxx_nml" is not found
    in "xxxx_xxxx.nml" (iostat=190).
                        ^^^^^^^^^^ ���줬��ħ.

���ξ��ˤ�, �ɤ߹��ޤ�ʤ� NAMELIST �����˶��Ԥ䥳���ȹԤ���������
�ʤɤ�, ���٥ץ����μ¹Ԥ�ư��Υ����å���ԤäƤ�������. 
���Ԥ䥳���ȹԤ򲿹Ԥ�­�����ȤǾ嵭�ξɾ��ϲ��򤵤��褦�Ǥ�.
(�����������Ǥ�). 

=end JA

=begin EN
== Fujitsu Fortran

=== NAMELIST loading error

A part of NAMELIST group name is not loaded in rare cases 
when a lot of NAMELIST group names is described in a file. 

Please check messages output by an execution program when information
of NAMELIST is not reflected in the execution program though it is
sure to be written in the file correctly.
The symptom described in the above-mentioned will appear if  
following messages are displayed.

 !*** WARNING [XxxxNmlRead] ***  NAMELIST group "xxxx_nml" is not found
    in "xxxx_xxxx.nml" (iostat=190).
                        ^^^^^^^^^^ This is the feature. 

In this case, please insert some null lines or comment lines in
front of the NAMELIST, and check operation with repeated
execution of the program.
The above-mentioned symptom seems to be evaded by some null lines or
comment lines. (The cause is uncertain).

=end EN


=begin HTML
<hr />
<small>
  $Id: compiler_note.rd,v 1.2 2014/07/08 07:47:12 sugiyama Exp $
</small>
=end HTML
