function cen_dis = get_center_dis(objs1,objs2)

    cen1 = [objs1(1)+1/2*objs1(3),objs1(2)+1/2*objs1(4)];
    cen2 = [objs2(1)+1/2*objs2(3),objs2(2)+1/2*objs2(4)];
    cen_dis = sqrt((cen1(1) - cen2(1))^2 + (cen1(2) - cen2(2))^2);

end