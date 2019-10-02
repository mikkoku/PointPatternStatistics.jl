window = (x=(-1.2,1.5), y=(5.4,6.2))
xy = [(0.171146, 6.18128), (0.825792, 6.1237), (0.900563, 6.13717),
    (-0.543602, 5.82919), (-0.629489, 5.40847), (0.864488, 5.59834),
    (0.0382095, 5.83211), (-0.573306, 6.16146), (-0.118548, 5.74562),
    (0.137276, 5.72065), (-0.387608, 5.45271), (0.709483, 5.68043),
    (-0.769612, 6.06163), (0.523729, 6.07474), (0.389607, 5.90159),
    (-0.559866, 6.02546), (1.31894, 5.84089), (-1.03956, 5.54941),
    (0.765956, 5.61121), (-0.246564, 5.52229), (0.860113, 5.52206),
    (0.576003, 6.12809), (-0.179952, 5.4613), (-1.14176, 6.01474),
    (0.689314, 5.58003), (1.49035, 6.00802), (0.363703, 5.40644),
    (-0.422678, 5.50465), (0.819573, 5.98695), (-0.08627, 6.05889)]
ref = [0.0, 4.3653322423953256e-5, 0.00017530848427371737, 0.0003954137396170454,
    0.0007023905710382694, 0.0010941330363524182, 0.0015666274222397858, 0.0021217332156657775,
    0.0027613239509604037, 0.003484672246775289, 0.0042913167206244385, 0.00517892470218595,
    0.00614875508238788, 0.007202760029804356, 0.008341160543495763, 0.009564503911010491,
    0.01087409358071667, 0.012269344698312201, 0.013750792125568512, 0.015320031173485349,
    0.01697293702895064, 0.018708352664438044, 0.0205289300602286, 0.0224348705029499,
    0.024428426110674306, 0.026509091137849383, 0.028677567773440216, 0.030931700826614672,
    0.03326768504771371, 0.03568943050206452, 0.0381917913940748, 0.04077314873895044,
    0.043419613944554625, 0.046125242533531674, 0.048900285364247065, 0.05174871917671808,
    0.05467257486607324, 0.05766438503293714, 0.06072045899357337, 0.06380869016929558,
    0.06693279226313487, 0.07010646532094111, 0.07332940871453009, 0.07659523992013828,
    0.07991266880424841, 0.08328269017007739, 0.08668307271402931, 0.09012249976969666,
    0.09360888067572704, 0.0971464403870087, 0.10073595750064934, 0.10436252694482595,
    0.10802437936502629, 0.11172157658825133, 0.11545979228092684, 0.11924955042305829,
    0.12309190992784913, 0.12698781546849924, 0.13093984596269448, 0.13494881691918492,
    0.13901412290209092, 0.14313832311172603, 0.14730187769406045, 0.15151181894848642,
    0.15577136564643812, 0.16008148010916545, 0.16444270743772182, 0.16885444111056125,
    0.17332000496449818, 0.17782367112728947, 0.18233519069991422, 0.18686459427444946,
    0.19141130085872982, 0.19598651251784016, 0.2005942384331224, 0.2052296419300934,
    0.20986440774184667, 0.21451658061013812, 0.2191921993089041, 0.22389011449437113,
    0.22859796856603665, 0.233321490175584, 0.23805854733009113, 0.24281739177584483,
    0.24760222733044834, 0.25241471041712726, 0.25725635303125516, 0.2621285633981171,
    0.2670314569532928, 0.27196404263915064, 0.2769194451978795, 0.28185621923197546,
    0.28677862873022897, 0.29170394795948806, 0.296639130861576, 0.3015888196792571,
    0.30655339227950207, 0.3115365099534235, 0.31653820593547866, 0.3215489326308548,
    0.3265528852257139, 0.33156035051241517, 0.3365763311491947, 0.34160407159663797,
    0.34664654645486126, 0.3517080323709918, 0.3567913155315914, 0.3618840544208992,
    0.3669715079865965, 0.37206580814734935, 0.37716060851673616, 0.38223965941515536,
    0.3873160020065506, 0.3923963874015237, 0.3974826481353848, 0.4025759181684342,
    0.4076840068459192, 0.4128089505607022, 0.41795047702792987, 0.4231107072721887,
    0.42829029452455913, 0.43349001884597405, 0.43870952913851835, 0.4439497426336315,
    0.44921155495659426, 0.4544956322902398, 0.45980184876528274, 0.4651322243064496,
    0.47048584836835483, 0.47585767445454086, 0.48120171305696635, 0.4865118725229576,
    0.4918088221886212, 0.4970974541522134, 0.5023820886054363, 0.5076654508620471,
    0.5129484758793467, 0.5182341208828642, 0.5235235399928748, 0.5288185540103508,
    0.5341199082311753, 0.5394267209091863, 0.5447222623341462, 0.550014491055092,
    0.5553086026585501, 0.5606068831140114, 0.5659102102335019, 0.5712204773380845,
    0.5765405320214001, 0.5818682021236018, 0.5872029706731531, 0.5925370793926497,
    0.597874518789654, 0.6032179903279467, 0.6085673763674821, 0.6139262675819204,
    0.6192715296607954, 0.6245868561243617, 0.6298842253115768, 0.635157089302571,
    0.6404167200697142, 0.6456516896333067, 0.6508498986099556, 0.6560242456690262,
    0.6611374893276698, 0.6662060215622143, 0.6712418835868752, 0.6762499414534344,
    0.6812314854721027, 0.6861885008878573, 0.6911227211197069, 0.696034368600978,
    0.700926829085495, 0.7058009550031068, 0.7106248478748473, 0.7153905520579944,
    0.7201155075415395, 0.7248072954285946, 0.7294672938522346, 0.7340980395065325,
    0.7386991352512237, 0.7432582943104953, 0.7477719163148111, 0.7522419705225791,
    0.7566727187886603, 0.7610648401995753, 0.7654218756658812, 0.7697421483832405,
    0.7739856495669473, 0.778184506387009, 0.7823436834703668, 0.7864598207658212,
    0.7905345343205491, 0.7945746459809095, 0.7985832393164524, 0.8025430530865227,
    0.8064171223231932, 0.8102328402692378, 0.8139989877398341, 0.8177104836416438,
    0.8213730377071582, 0.8249878569708711, 0.8285574417979189, 0.832085355982109,
    0.8355720068056933, 0.8390199341387028, 0.8424303607083435, 0.8458029036641665,
    0.8491392263048995, 0.8524396287622842, 0.8557033083682104, 0.8589410390162966,
    0.8621536821033504, 0.865340601334038, 0.8685036668628742, 0.8716416392240589,
    0.8747529132788527, 0.8778165889508787, 0.8807810162479123, 0.8836276129663937,
    0.8864137653591784, 0.8891544555771238, 0.8918512759569056, 0.8945039460548612,
    0.8971066495898357, 0.8996630797540144, 0.9021779793636634, 0.9046640238455361,
    0.9071234433997776, 0.9095589221085179, 0.9119701262004352, 0.9143569373502427,
    0.9167205830294206, 0.9190624170457723, 0.9213817664563575, 0.9236795960482577,
    0.9259557174574908, 0.9282113217176574, 0.9304456554011953, 0.9326587737902291,
    0.9348529649639792, 0.9370262524052397, 0.9391817122886855, 0.9413185588741494,
    0.9434366007076981, 0.9455364521281238, 0.9476183397424743, 0.94968309837709,
    0.9517297171437048, 0.9537165413323675, 0.9556599963451246]
r = 0:0.001:0.25
est = Fest(xy, window, r)
@test maximum(abs, (est - ref)) < 0.01
est = Fest(xy, window, r, 1000)
@test maximum(abs, (est - ref)) < 0.001
