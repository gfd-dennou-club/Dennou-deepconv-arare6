=begin JA

= ��ʿ��ľ 2 �����η׻���Ԥ��ˤ�

# * �����̰�ϯ (sugiyama)
#   * $Id: calc_2d.rd,v 1.1 2014/03/01 20:11:21 sugiyama Exp $


��ʸ��ǤϿ�ʿ��ľ 2 �����η׻���¹Ԥ��뤿�����ˡ�򼨤�. 


== ����ե�����ν���

����ե�����(NAMELIST �ե�����) �� &gridset_nml, &axesset_nml ���ѹ���
������Ǥ���. (({ydim = 1})) �Ȥ��뤳�Ȥ�, deepconv/arare5 �Ͽ�ʿ��ľ
2 �����Υ�ǥ�Ȥ���ư���. 

    &gridset_nml
      xdim  = 50                 ! X �����������
      ydim  = 1                  ! Y �����������
      zdim  = 50                 ! Z �����������
      NCMAX = 1                  ! ���ؼ�ο�
    /

    &axesset_nml
      Xmax  = 1.0d3             ! X ��ɸ�ν���
      Zmax  = 6.5d2             ! Z ��ɸ�ν���
    /

���ΤȤ�, Ymin, Ymax �ϻ��ꤹ��ɬ�פ�̵��. 


=end JA


=begin HTML
<hr />
<small>
  $Id: calc_2d.rd,v 1.1 2014/03/01 20:11:21 sugiyama Exp $
</small>
=end HTML

