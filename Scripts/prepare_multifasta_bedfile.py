
# the overall point of this script is to make a bed format file containing the 
# positions to be masked from the genome the output format is 
# contig_name	start_position	end_position
# this script is run, and the output piped to a file e.g. 2022.09.06.tb_genomes_to_be_masked.bed
# then bedtools is run
# bedtools maskfasta -fi my_sequences.fasta -bed 2022.09.06.tb_genomes_to_be_masked.bed -fo my_sequences.masked.fasta

import sys

# samples is a list of all the contig names that need to be masked
samples = sys.argv
del samples[0]
print(samples)

# positions_to_mask is a list of tuples, where each tuple is a range within the reference
# genome that needs to be masked
positions_to_mask = [(23173, 23273), (33582, 33794), (79507, 79551), (80236, 80550), (103713, 105215), (105324, 106715), (131382, 132872), (149533, 150996), (154073, 154125), (154126, 154178), (154179, 154231), (177543, 179309), (179319, 180896), (187433, 188839), (188931, 190439), (206812, 206850), (206869, 206907), (272855, 272955), (307877, 309547), (333437, 336310), (336560, 339073), (339364, 340974), (349624, 349932), (349935, 351476), (361334, 363109), (366150, 372764), (372820, 375711), (399535, 400050), (424269, 424694), (424777, 434679), (466672, 467406), (467459, 468001), (530751, 532214), (543174, 544730), (547488, 547517), (580578, 580654), (580655, 580731), (580732, 580808), (616828, 616878), (622793, 624577), (642754, 642811), (671996, 675916), (701247, 701369), (703912, 703985), (706790, 706863), (709425, 709548), (709585, 709663), (711624, 711702), (795467, 795518), (802429, 802477), (812835, 812921), (812922, 812975), (832534, 832848), (832981, 833508), (835701, 838052), (838451, 840856), (846159, 847913), (848103, 850040), (850342, 850527), (863155, 863255), (889017, 889020), (889021, 889048), (889072, 889398), (889395, 890333), (890348, 890375), (890376, 890379), (890388, 891482), (921575, 921865), (924951, 925364), (925361, 927610), (927837, 930485), (947312, 947644), (960173, 960225), (960226, 960278), (960279, 960333), (968424, 970244), (976872, 978203), (1020058, 1021329), (1021344, 1021643), (1025458, 1025472), (1025497, 1026816), (1026879, 1026893), (1027061, 1027076), (1027104, 1027685), (1027685, 1029337), (1029345, 1029360), (1090373, 1093144), (1093361, 1094356), (1095078, 1096451), (1158918, 1159307), (1159375, 1160061), (1160095, 1160433), (1161297, 1162472), (1162549, 1163376), (1164572, 1165435), (1164572, 1164589), (1165092, 1165499), (1165532, 1165549), (1169423, 1170670), (1176928, 1177242), (1177239, 1177373), (1179345, 1179395), (1188421, 1190424), (1190757, 1192148), (1211560, 1213863), (1214513, 1214947), (1214769, 1215131), (1216469, 1219030), (1251621, 1252945), (1262272, 1264128), (1276296, 1277643), (1277843, 1277846), (1277847, 1277863), (1277893, 1278300), (1278800, 1278816), (1278817, 1278820), (1298764, 1299804), (1299822, 1300124), (1301755, 1302681), (1305495, 1305556), (1305557, 1305618), (1305619, 1305661), (1339003, 1339302), (1339349, 1340524), (1341358, 1342605), (1357293, 1357625), (1384989, 1386677), (1456585, 1456627), (1457453, 1457504), (1457505, 1457557), (1468143, 1468161), (1468171, 1469505), (1469633, 1469651), (1488154, 1489965), (1507531, 1507581), (1532443, 1533633), (1541949, 1541951), (1541994, 1542980), (1542929, 1543255), (1543307, 1543309), (1561464, 1561772), (1561769, 1563388), (1572127, 1573857), (1606386, 1607972), (1612558, 1612578), (1612579, 1612599), (1612600, 1612620), (1612621, 1612641), (1612642, 1612662), (1618209, 1619684), (1625366, 1625418), (1630638, 1634627), (1633531, 1634790), (1636004, 1638229), (1637133, 1638392), (1644261, 1644313), (1644314, 1644364), (1655609, 1656721), (1751297, 1753333), (1779266, 1779277), (1779314, 1779724), (1779930, 1780241), (1779959, 1780047), (1780199, 1780699), (1780485, 1780573), (1780643, 1782064), (1782072, 1782584), (1782758, 1783228), (1783309, 1783623), (1783620, 1783892), (1783906, 1784301), (1784497, 1785912), (1785912, 1786310), (1786307, 1786528), (1786584, 1787099), (1787096, 1788505), (1788162, 1789163), (1788514, 1789811), (1788514, 1788525), (1830074, 1830125), (1855764, 1856696), (1862347, 1865382), (1907460, 1907515), (1907516, 1907571), (1927218, 1928589), (1931497, 1932654), (1932694, 1933878), (1938093, 1938145), (1944756, 1944808), (1981614, 1984775), (1982887, 1982964), (1982965, 1983042), (1983043, 1983120), (1983121, 1983198), (1983199, 1983276), (1983277, 1983354), (1987703, 1987730), (1987745, 1988629), (1988680, 1989006), (1989030, 1989057), (1989833, 1992577), (1996101, 1996128), (1996152, 1996478), (1996529, 1997413), (1997428, 1997455), (1998584, 1998597), (1999142, 1999357), (1999800, 1999813), (2000614, 2002470), (2025301, 2026398), (2026477, 2026776), (2026790, 2027971), (2028425, 2029477), (2029904, 2030203), (2039453, 2041420), (2042001, 2043272), (2043384, 2044775), (2044923, 2046842), (2048072, 2048371), (2048398, 2049597), (2049921, 2051150), (2051282, 2052688), (2059441, 2059498), (2059518, 2059575), (2061178, 2062674), (2087971, 2089518), (2162932, 2167311), (2163323, 2163392), (2163393, 2163461), (2163462, 2163530), (2163741, 2163809), (2163810, 2163878), (2163879, 2163947), (2163948, 2164016), (2164017, 2164085), (2167649, 2170612), (2195989, 2197350), (2226244, 2227920), (2260665, 2261144), (2261098, 2261688), (2330147, 2330225), (2343027, 2343332), (2356729, 2358033), (2365414, 2365441), (2365465, 2365791), (2365788, 2366726), (2366741, 2366768), (2367359, 2367655), (2367711, 2368442), (2372437, 2372492), (2372494, 2372549), (2381071, 2382492), (2387202, 2387972), (2423240, 2424838), (2430117, 2430144), (2430159, 2431199), (2431094, 2431420), (2431444, 2431471), (2439282, 2439947), (2458392, 2458449), (2493801, 2493818), (2522173, 2522230), (2523184, 2523236), (2531898, 2531950), (2531951, 2532003), (2532004, 2532056), (2532057, 2532109), (2532110, 2532162), (2532163, 2532212), (2550011, 2550013), (2550014, 2550041), (2550065, 2550391), (2550388, 2551326), (2551341, 2551368), (2551369, 2551371), (2581843, 2582298), (2600731, 2601879), (2617667, 2618908), (2632923, 2634098), (2634528, 2635592), (2635577, 2635604), (2635628, 2635954), (2635951, 2636889), (2636904, 2636931), (2637688, 2639535), (2651753, 2651938), (2687128, 2687179), (2687180, 2687257), (2692799, 2693884), (2706017, 2706736), (2716315, 2716391), (2720644, 2720656), (2720776, 2721777), (2721844, 2721856), (2727336, 2727920), (2727967, 2728266), (2762762, 2763061), (2763397, 2763696), (2784614, 2784642), (2784657, 2785697), (2785592, 2785918), (2785942, 2785970), (2795301, 2797385), (2800671, 2800918), (2801254, 2806236), (2806368, 2806625), (2828556, 2829803), (2835785, 2837263), (2921551, 2923182), (2935046, 2936788), (2943600, 2944985), (2960105, 2962441), (2970551, 2971549), (2972106, 2972108), (2972109, 2972136), (2972160, 2972486), (2972435, 2973421), (2973436, 2973463), (2973464, 2973466), (2973795, 2975234), (2975242, 2975775), (2975928, 2976554), (2976586, 2976909), (2976989, 2977234), (2977231, 2978658), (2978660, 2979052), (2979049, 2979309), (2979326, 2979688), (2979691, 2980818), (2983019, 2983033), (2983071, 2983874), (2996003, 2996053), (2996054, 2996104), (2996105, 2996155), (3007063, 3007115), (3007116, 3007168), (3007169, 3007221), (3013612, 3013687), (3053914, 3055491), (3073055, 3073112), (3076894, 3078078), (3078158, 3078985), (3079309, 3080457), (3100202, 3101581), (3115741, 3116142), (3116818, 3118227), (3119185, 3123576), (3119185, 3119220), (3119259, 3119294), (3119335, 3119370), (3119411, 3119446), (3119484, 3119519), (3119556, 3119591), (3119627, 3119662), (3119701, 3119736), (3119777, 3119812), (3119848, 3119883), (3119921, 3119956), (3119995, 3120030), (3120068, 3120103), (3120141, 3120176), (3120213, 3120248), (3120285, 3120320), (3120359, 3120394), (3120433, 3120468), (3120504, 3120523), (3120566, 3121552), (3121501, 3121827), (3121882, 3121897), (3121938, 3121973), (3122013, 3122048), (3122086, 3122121), (3122158, 3122193), (3122230, 3122265), (3122303, 3122338), (3122375, 3122410), (3122436, 3122471), (3122513, 3122548), (3122585, 3122620), (3122661, 3122696), (3122738, 3122773), (3122811, 3122846), (3122882, 3122917), (3122955, 3122990), (3123029, 3123064), (3123102, 3123137), (3123173, 3123208), (3123248, 3123283), (3123318, 3123353), (3123390, 3123425), (3123467, 3123502), (3123541, 3123576), (3155874, 3155927), (3155928, 3155981), (3155982, 3156035), (3156036, 3156089), (3160522, 3160583), (3162268, 3164115), (3171468, 3171518), (3171522, 3171572), (3171576, 3171616), (3181794, 3181836), (3192202, 3192254), (3192255, 3192307), (3192308, 3192360), (3194166, 3195548), (3200794, 3202020), (3288464, 3289705), (3289705, 3290235), (3289790, 3290506), (3313283, 3313672), (3318835, 3318889), (3319468, 3319568), (3319569, 3319666), (3333768, 3333773), (3333785, 3335164), (3335787, 3335792), (3376939, 3378243), (3378329, 3378415), (3379376, 3380452), (3380440, 3380682), (3380679, 3380993), (3381351, 3381365), (3381375, 3382622), (3382660, 3382674), (3465778, 3467091), (3481399, 3481413), (3481451, 3482698), (3482708, 3482722), (3490476, 3491651), (3501334, 3501732), (3501794, 3502936), (3510088, 3511317), (3527391, 3529163), (3551227, 3551229), (3551230, 3551257), (3551281, 3551607), (3551604, 3552542), (3552557, 3552584), (3552585, 3552587), (3552710, 3552712), (3552713, 3552740), (3552764, 3553090), (3553087, 3554025), (3554040, 3554067), (3554068, 3554070), (3557311, 3558345), (3591493, 3591569), (3626614, 3626666), (3658658, 3658715), (3704895, 3705004), (3710382, 3710409), (3710433, 3710759), (3710756, 3711694), (3711709, 3711736), (3711749, 3713461), (3729364, 3736935), (3736984, 3738438), (3738158, 3742774), (3743198, 3743404), (3743402, 3743510), (3743508, 3743605), (3743711, 3753184), (3753765, 3754256), (3754293, 3755237), (3755952, 3767102), (3769514, 3769720), (3769754, 3769862), (3770994, 3771091), (3778568, 3780334), (3795058, 3795085), (3795100, 3796086), (3796035, 3796361), (3796385, 3796412), (3799987, 3800011), (3800092, 3800796), (3800786, 3801463), (3801530, 3801554), (3801653, 3803848), (3842239, 3842769), (3843036, 3843734), (3843885, 3844640), (3844738, 3845970), (3847165, 3847701), (3847642, 3848805), (3849294, 3850139), (3883550, 3884921), (3890779, 3890806), (3890830, 3891156), (3891051, 3892091), (3892106, 3892133), (3894093, 3894389), (3894426, 3895607), (3926569, 3930714), (3931005, 3936710), (3939617, 3941761), (3941724, 3944963), (3945098, 3945597), (3945794, 3950263), (3950830, 3951329), (3969343, 3970563), (3970705, 3972453), (3978059, 3979498), (3991568, 3991625), (3997980, 3999638), (4031404, 4033158), (4036731, 4038050), (4052949, 4052966), (4052971, 4052994), (4052995, 4053105), (4053004, 4053021), (4053106, 4053216), (4053217, 4053327), (4053328, 4053438), (4053439, 4053549), (4060648, 4061889), (4061899, 4062198), (4075615, 4075630), (4075752, 4076099), (4076484, 4076984), (4076984, 4077730), (4077735, 4077750), (4078506, 4078518), (4078520, 4079749), (4079786, 4079798), (4091233, 4091517), (4093632, 4093946), (4093940, 4094527), (4134601, 4134725), (4189285, 4190232), (4190284, 4190517), (4196171, 4196506), (4198874, 4199089), (4215881, 4216063), (4252993, 4254327), (4276571, 4278085), (4301563, 4302789), (4302786, 4303397), (4318775, 4319266), (4348721, 4348773), (4348774, 4348826), (4350745, 4351044), (4351075, 4352181), (4353280, 4353330), (4353331, 4353381), (4353382, 4353432), (4374484, 4375683), (4375762, 4375995)]

# simple for loop to print the file
for sample in samples:
	for window in positions_to_mask:
		print('\t'.join([sample, str(window[0]), str(window[1])])) 
