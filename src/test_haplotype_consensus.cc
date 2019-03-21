#include <iostream>
#include <sglib/workspace/WorkSpace.hpp>
#include <sglib/processors/LinkageUntangler.hpp>
#include <sglib/processors/HaplotypeConsensus.hpp>



int main(int argc, char **argv) {
    WorkSpace ws;

    ws.sg.load_from_gfa("initial_graph.gfa.gfa");

    ws.long_read_datastores.emplace_back();
    ws.long_read_mappers.emplace_back(ws.sg,ws.long_read_datastores[0]);
    ws.long_read_mappers[0].read_filtered_mappings("fm_10K3.bsgfrm");
    ws.long_read_mappers[0].update_indexes();
    ws.long_read_mappers[0].improve_filtered_mappings();
    ws.long_read_mappers[0].read_paths.resize(ws.long_read_mappers[0].filtered_read_mappings.size());

    ws.linked_read_datastores.emplace_back();
    ws.linked_read_mappers.emplace_back(ws.sg,ws.linked_read_datastores[0],ws.uniqueKmerIndex,ws.unique63merIndex);
    ws.linked_read_mappers[0].read_tag_neighbours("tag_neighbours_1000_0.03.data");



    auto u = LinkageUntangler(ws);
    LinkageDiGraph imldg(ws.sg);

    LinkageDiGraph ildg1(ws.sg);

    LinkageDiGraph emptyLDG(ws.sg);

    std::vector<std::vector<sgNodeID_t>> lines {
            {
                -2186685, 34190, 2129381, -2021024, -2144860, -1922595, 2131337, -2090596, 1844227, 1980845, -1857288,
                -2165057, -2119315, -2094320, 1953923, 2163069, 1493149, -2009038, -2086713, 1832192, -2200423, -2107299,
                -2128970, 2011369, 2177332, -2072383, 2151835, 2099154, 2063234, 2167743, 2184527, 2042763, 2201839,
                -2183086, 1981004, 1861648, -2008873, -1940556, -2141561, 2149023, -409488, -2139718, 1961464, -2162023,
                -2028484, 2119349, 2172266, -2175361, -1849820, -2087089, 1926194, -2155299, 2134678, 1909169, -1944934,
                -2155162, -2167896, 1525457, 2173980, -2115298, -2178200, 1952466, 2185487, 2104071, -2107388, 1954525, 51,
                2171773
            },//0
            {
                1979694, 2021015, 1997321, 2146399, -1984150, 1992936, -2047239, 2086799, 86, -1895332, -2184470,
                -1992989, 2083318, -1953893, -1249150, -2123007, -2168737, 1969892, -1912396, 2187689, 2077498, 2139013,
                1903994, 1027970, -1573256, 1986997, 1961485, 2067813, 2204521, 1917185, 1934438, -2003786, 2133317,
                -1874989, 1970631, -1918803, 1996416, 1940268, -2140412, 2137187, -2128258, -2181978, -1991742, -1954614
            }, //1
            {
                2022498, 2079919, -2152872, -2011878, -2106342, -2091070, -2104479, -2030108, -2020564, -1989530, 2096213,
                -2044376, 2051186, 1995979, 2178258, -1932239, -1958795, -1951544, -2127385, -1999656, -2099350, -2166102,
                2148854, -2156141, -1960738, 2028430, -2128106, -2194262, -1044039, 1977974, -639499, -1967399, 1957618,
                -2065561, 2184426, 2157223, -2005145, -2189711, -947321, 1952037, -1855437, -1903414, 1935956, 1937911,
                -2091221, 2075440, -1841319, -1990031, 2041975, 2028076, -1980874, -2141221, 2028321, -1870433, -2148960,
                2067660, 2203770, -798201, -2147695, 2094388, 2026255, -1945224, 1973784, 1962854, 2112511, -1978612,
                2123990, -2100682, 2180605, -1878404, 2097627, -1992542, 1877445, 2175970, 2063832, -2012493, 111,
                2028523, -2022515, -969771, -2158310, -2075139, -2115094, -2042908, -2019343, 2169343, -1493403, 2100561,
                -2131071, 2103558, 2150750, 2161302, 2082796, -2075155, -1962993, -1981423
            },//2
            {
                2060639, -2114550, 147279, -2188835, 1722144, -2083054, 1978466, 2113986, -1836318, 1944500, -2186452,
                2125654, 2075615, 294, -1947153, -1856709, 1925103, -1899060, -1952198
            },//3
            {
                1585926, -384299, 2054116, 2164877, 2150711, -2200232, -2149847, 2167685, 1987131, -1832700, -2022418,
                -2122499, -2049652, 2017363, 1333314, 2000985, 1989440, 1903445, 1827912, -2178998, 1823244, 2101648,
                2144079, 2062282, 1901200, -2076346, -2108199, -2130769, -1926785, 1427356, -2168248, -2099629, -2043192,
                -2053159, 1919288, -1990118, 1883913, -2111505, -2177075, -1882937, -2124318, 2026877, 2027412, 2178224,
                -2102755, -1926959, -1009046, 1238615, 1946749, -1495749, 20446, -2205975, 1918209, -2172625, -2074821,
                1880793, 2162542, -2011787, 1923733, 2069456, 164190, 2203766, -2037673, 1895690, -1973907, 2076080,
                -2149630, -1891353, -2148329, -2084819, -1912309, -2182482, 2081188, -488734, 345047, 1971311, 1857236,
                1621705, -2160971, 1871814, -1948336, 1583605, 2181630, -1912845, -2101518, -2170008, 1969092, -2065066,
                -1978874, -2153263, 2026723, 1525519, -2139981, -1492077, -2183999, 2047550, -2068413, -1866561, 1882043,
                1935845, 1879835, 2200413, -2198830, 2014387, 1933113, -1767018, -2084880, -1832050, 1975399, -2169508,
                -2002096, 1964308, 373, 2053498, 2118813, -1957322, 1907004, 1966548, -1816873, -1856357, 2021730, 2124590,
                1869010, -2176342, -1549856, 1989264, -1999194, -2030545, -2157622, 1994474, 2166783, 1104415, -1993398,
                2145325, -2156861, -2196474, -2097490, 1942705, 2155662, -2044926, -2075691, -2094695, 2092614, 2088318,
                -2153561, 2092162, 2165402, -2082126, 1859367, 2058113, -1982781, 2176075, 2205047, -1923788, -2116439,
                2121988, -1817817, -2133958, 2177341, -1916902, -2157732, 1947368, -607618, 1995284, -2113979, 2153790,
                1913044, -1856524, 2031463, 1853997, 1902020, -1574537, -2076250, -2039056, -2127909, 2173617, -1833508,
                -1972317, -1924989, 2142906, -1929231, 2133967, 2005287, 1948670, -2021470, -2186568, 2007125, 2173840,
                -1977039, 2201386, 1904754, -2068385, 2175205, 694677, 2099781, 1878742, 1600393, -2130127, -2125479,
                -2168761, -652196, -2131200, -2151563, 2188229, -1956184, -2061501, 2098617, 2196217, 2160335, -2190425,
                2059534, 2194073, 2113214, 2068721, -1822855, -106282, -1426140, 2187448, -1961695, 2162755, -2134324,
                -2100339, 2128239, -1491135, -1556521, 2131065, 2014578, -1872173, 2011937, -2132991, -2126333, 1985354,
                -2144096, 208080, 1941509, 2201573, 2160037, -2098589, -1982320, 2115484, 2101775, 2147468, 1961050, 1743672,
                2079190, 2127055, 2090861, 2171600, -1919235, 2069022, 1933434, -2158607, -2033506, -1778663, -2041588,
                2017897, 1939102, -2099676, -1875348, -2084104, 1890615, 2121190, -2087887, 2117917, 2060947, 1951956,
                -1973929, -1975250, -1855895, 1927624, -1993221, 2081973, 2144806, 2125870, -2051715, 1903847, 1686893,
                2132566, 1440182, -2082747, 2173725, -1978683, -2131865, -1936959, -2095028, 2043639, -2159047, -2045654,
                -2161803, -2104528, 2190156, -2186144, 519929, 2061744, -2075721, -2013763
            },//4
            {
                -1450272, -1824780, 1970057, -1899053, -2009584, 1914415, -2078342, -1999672, -2025784, -1980637,
                2085743, 2011554, 2087262, -2162500, 1483644, 2122529, 1427848, -2119173, -2165206, -1983772, -1862988,
                2127134, 1721334, 1898517, 413, -92314, 2042078
            },//5
            {
                2067848, -2042005, 2148824, -838843, -2051938, 1978141, -1582629, -2206801, 1895590, -1899488, 1614357,
                1886556, -2006019, 1311774, 1564207, -1554425, -1844270, 1895908, -1871789, 2054426, 1900587, -1884048,
                -1845725, 1839364, 1947756, 2162895, -1566316, -1926072, 1997554, 2191458, 1894707, -2124426, 1866024,
                1995762, -2018722, -2194755, -2163083, 525, -2098453, 1569473, -1475223, -1848954, -1531142, 1627303,
                2112135, -2019023, 345737, -2077338, 1987198, 2182708, -1907192, -1937935, -2068758, 1833431, 2196763,
                -1906674, -2156181, 1988144, 2014874, -1976721, 1898232, 1841407, -2134123, -2037945, -1523915, 2058427,
                1899972, -1864953
            },//6
            {
                -1555186, -2195393, -1934963, -1008365, 1883831, 2113758, -1932103, -1906469, 2006678, -2066197, -1618673,
                1477450, 2116373, -1722691, -1952172, -2124864, -1258908, -2112436, -2049642, 1925073, 1490105, -1402027,
                -1967375, -2038420, -2197969, -1820011, -1980043, -734557, -2234336, -1514125, 1873910, 1498246, 1593445,
                -272779, 2076441, -1241106, -1888122, -2022145, -1853566, -1985309, 2189988, 1886587, 2100118, 1520306,
                2200782, 1628408, 454912, 1951102, -2135661, 1879705, -2095109, -1999446, 1437955, -1516182, -2019570,
                1933813, -2206802, -2194529, 2064972, -2159580, -1940075, 2204970, 1512320, 615807, -1569987, 1876096,
                -2133923, 2188427, 1481337, 1723056, -2007957, -1929395, -1918651, 1621945, 1570175, 2157676, 1508726,
                -1456291, -2031231, -2128283, -1433848, -1573796, 1547174, 1464005, 2068325, -2155303, -1512644, -2096651,
                -2115364, 1486066, -1896927, 2231914, 1985684, -2004848, -2081239, 2203644, 1910865, -1883048, 1840229,
                1711364, 2170587, 2097079, 1888036, 2039607, -1519708, -2105054, 1934433, -2145996, -1996615, -2075138,
                2206000, 1870909, -2139467, 1982962, 1798656, 1864049, 2064854, 1555767, 2204941, -1833833, 1545006,
                1935582, -1549748, 2018288, -267919, 2193653, -1560960, 2092820, 1887660, -2193200, 1463167, 1669105,
                -1927322, -1902740, -1504948, 2163818, 2192801, -2036325, -2166790, -1571618, -2139290, 2055109, -1998703,
                -1945915, 2103995, 2194146, -1953792, 1552174, 515247, -1973592, -2006713, -2105709, 1900650, 1567926,
                2204430, -2234166, -1878551, -2011125, -1836975, 2072015, 1565153, 2214353, -1544976, -2045492, -413674,
                -1280336, 2161853, -1981453, -1723761, 1987491, -1933522, 1992185, -1555162, 1521258, 1499193, -1883692,
                -2198477, 2152475, 2101062, -1524804, -2026651, 2039506, 1855715, -1412616, 2189618, 2009459, -1828820,
                -1457964, 1940908, 1988743, -1872056, 1567524, 2013726, -2006444, -1500060, -1534300, 2039109, 2085378,
                -2178297, -1525250, 1532937, 1955627, -834475, 1992111, 1851640, 548, 1922049, 1674348, -1987101, -1888931,
                1582104, 2024460, 1927669, 1469646, 1914522, -1582727, -1953898, -1845353, -1837323, -2040969, -1786569,
                2038041, 1934725, 2051573, 1478714, 2153196, -2006222, 2048306, -1913591, -2037076, 2183949, 1896765,
                1995119, -2188151, 1565462, 1862346, 2000341, 1859323, -1926409, 1419438, 1877835, 2026154, -2185743,
                -2015960, -1486051, 1591297, 1893736, -1446581, -1566030, -1536374, -1916936, 1561738, 2042343, -1952933,
                2144407, 2066761, 2137352, 2171259, 1572251, -1566506, -1992360, -1867617, -1910944, 1989987, -1582665,
                -1965303, 1959190, -1566852, -1856967, -2132382, 2142436, -1917507, 1561275, 2162740, -1567595, 2079604,
                1972795, 1628200, 2095956, 1840669, -1567496, 1872322, -1998638, -1872686, 1829446, -2130472, -2036284,
                -2008055, 2045375, -521017, 1828337, 1843353, -1510225, -2109610, -1475900, 1995712, -2104112, -26029,
                -1895489, 2054546, -2097714, 2189764, 439505, -1977712, -1799046, 2148665, 2186011, 1579549, 1718756,
                1933059, -1569834, -2222610, 828133, -2084159, -1845917, -1991671, -2201008, -1534644, -1581771, 2081811,
                1870730, 1432607, 1576705, 1550748, -2148205, 1395084, -1912313, -2109402, -2005088, 1457571, 2155174,
                1986404, 2114054, 1569543, 1862995, -1471645, -1958942, -1829002, -1722526, -1965000, -2182839, -1873587,
                -2217082, 1847310, -1870383, -2120017, 2006462, -1823108, -1875400, -655983, 1991369, 1429452, 1931829,
                1356451, -995423, -2057513, 1427533, 2076838, -1911780, 1575875, -1637922, -1949322, -1877334, -1551083,
                1896063, 2190476, 1483406, -2049400, 1856418, -2082255, -2048393, 1960137, 1965244, -2202822, 1825330,
                1887386, 1988778, -2034753, -2088178, 1990670, 2129337, 1676440, 2146140, 2078934, 1441868, -1871681,
                2135130, 1757772, 1571168, 1721563, 1825043, -2072947, -2048821, 1517004, 1866696, -1950243, -1992266,
                1449907, -1923295, 1944832, 1615751, 2183804, -1919051, 1561038, 2074749, 2016967, -236291, 1542851,
                -1428848, -2152558, 1834020, 1873598, -2133574, 1948353, 2004495, 1577888, -2148065, 1879287, -1997239,
                1562490, 2069348, 1923637, 1495879, 2151531, 2196868, 1983535, -1556420, -1910413, 1526486, 2142470,
                -1996008, 1863216, 2003932, 2008361, -2195048, -2019792, -2072043, -1869059, -1455230, -2019913, -1534362,
                2044977, 999882, 2014169, 2205082, -1570393, -1676183, 1525065, 1958500, 1549859, -2136136, -1244809,
                1521957, 2019997, 1473861, 1903946, 1844261, -1509958, -2035344, 1960281, -2134216, -2160611, -1610396,
                1563558, 2156229, -972471, -757065, -1944703, 2121778, -979530, 2144580, -1514715, -1321772, -499648,
                1990289, -1957325, -1913304, 2048976, -1313230, -2134609, -1932145, 1839825, -1884315, -1966602, 1559743,
                -1558234, -2047642, -1634775, 1011310, -1975619, -1867562, 2165599, 2179828, -1913860
            },//7
            {
                2073624, -2096217, -2197615, 2062892, 1998675, -1995731, 311460, 2110391, 1266363, 2016882, 1955226,
                1950089, 2129963, 2039068, 2076414, -1777994, 2067312, -106233, -2185668, -1999365, -2044847, 2205555,
                -2182344, 2091589, -2054111, 2172375, 2016861, -2135693, 2142234, -1838170, 2084936, 2135422, -2144078,
                -1850408, -2103569, -2065909, -2178146, -2119064, 2177423, -2091255, -1966566, -2019678, -2003552,
                2010015, 2118701, 2119551, 2053470, -2094859, -2133926, -2166615, -1428889, -2178118, 2195554, 1995503,
                2095405, -2180131, -1852103, -2144861, -1940345, 1927099, -1836541, 1900074, -2179904, 2174611, -1918399,
                2039467, -2056134, 2007460, 2066531, 2115770, 2113583, 2135125, 2018962, 2140891, 1958510, 620997,
                -2201813, -2093710, -1913882, 2199386, -1958376, -2135570, -1897466, 1923483, 1944459, 2125612, -2033819,
                -2151569, 2182831, 1904089, -2128763, 1943769, 2130201, -2022859, -1868822, -2196305, -1852197, 2162605,
                -1966950, -469881, 1956756, -2090958, -2167439, 1966784, 2001950, -1989342, 2008786, 2142034, -2026867,
                1978099, -2060271, -2155886, 1950725, -2062101, -1829473, 1832830, -1931902, 1997444, -1916188, -2074883,
                -2099968, 2175718, -2057285, 1911313, 2194543, -1920356, -1931393, 1910123, -2048753, 2204739, -2041838,
                2103323, 2142865, -817591, -1912717, -2093204, -2197620, 2061146, -2170020, 2090074, -1916559, -358681,
                -2190707, -191307, 2175065, -2099566, 2087016, -2146436, 2034249, -2117893, -1984055, -2024250, 2147317,
                2174440, 1917949, 2181573, 2000093, 2117607, 2137993, 2049533, 1955239, 2113539, -2013662, -2062641,
                -2103184, 2202644, 2126155, 2169430, -2085358, -2118153, -2050535, 1926349, 1978231, -1950222, -2170586,
                -1889000, 1485502, 1905832, -2002041, -2204410, 2079775, 1930762, 1823546, 2202144, 2195444, -2179989,
                2012705, -2090726, -2077347, 2196940, 1547582, 1984094, 2124661, 2085033, -2117309, -1964259, 2141049,
                -2195445, 2112333, 2068800, 2079071, 2015314, -2065124, 1939614, -2068288, 2168446, 2156427, -2029508,
                1922832, -1906539, 2206758, 2129481, -2156984, 2087488, 2051315, 2009419, 733, -2184351, -1867243, 2170926,
                2112419, -2124630, 1894187, 1994734, 2058237, -1922439, 2052091
            },//8
            {
                1862135, 2201320, 2072118, -2147622, 2206732, 2234256, 1506034, 1580387, -1350592, -1700304, -1765492,
                -2157511, -1678881, -2185483, -2017663, 1956728, -2136258, 1763547, -2105171, -1581968, 2170844, 1575809,
                -1402462, -2177747, -937436, -2082962, -2104736, -2079057, -1544259, -1461143, 2224119, -1951319, -1918026,
                -1508543, 1982731, 2071794, 2141331, -1564984, 1872605, 1848616, -1832318, -2073353, -1998342, -1505742,
                2016964, 902, 2149046, -2133101, -1941496, -2149061, 1835687, -1918371, -302208, -540800, 1890036, 2164455,
                -1555465, -1588502, 1978074, 1881379, -1626580, -1555614, -2078499, 1452271, 1954634, -1874202, -1560582,
                2002384, 2049148, 2071098, 2035128, -2074145, -1494417, -1883146, 2068343, 1953917, 1610733, -2071779, 1547526,
                1880866, -2011170, -1886612, 2126867, -2070878, 2167267, -2163348, 2032122, 1887822, -2145761, -1871517,
                -2138663, 2070464, 2135992, -2082809, 1913067, 2130971, 1903930, -1882056, -1457150, -2032989, -1856219,
                1644097, -1938775, 1883938, -2074937, -2104277, 2147841, -1882082, -2074182, -2156930, 2107550, -2020799,
                -1982076, -2179693, -2166917, 2165377, -2081184, -1901606, -1997855, -2025565, -2120428, -2040280, -1933072,
                -1920563, 2025717, -1916166, 1865257, -2025466, 2033225, 1922885, -530474, 2034307, -1854495, 1963792,
                1919094, -2046268, -2027116, -2042660, 2000442, 2187490, 2148867, 2018094, -2163803, 2119346, 2166224,
                2098150, 2066170, -1997918, -1896911, 1914032, -2196066, -2204603, 2070855, 2146476, 2051519, -2092709,
                -1917705, -1970940, -2063308, 1869784, -1927980, -2099678, 1448565, 1919305, -1963133, -2133943,
                2156251, 2133194, 2118484, -2117736, 1939489, 2085623, 2106522, 1511885, -1614655, -1602536, 1465978,
                -1788139, 1272211, 2117954, -2114869, -1905661, 2061159, 1930702, 1933815, -2089012, 2002336, 2100476,
                2084422, 2111559, 1900873, 2157814, -1871249, -1536512, 1594191, -2166672, -1893199, -1908205, -1995157,
                -2181702, 2165333, 2201072, -2166209, 2055033, 2129581, -2199521, -2070931, -2088984, 2078388, -2167419,
                -1966677, -2059989, -2197630, 2162418, 2026243, 1936489, 1656561, 2182492, -1950667, 1924456, -2094192,
                -1857758, 2191365, 2166232, -2149988, -2152432, 2205129, -2121550, 1998628, -2037012, 2064674, -184425,
                2182301, -2078777, 2093908, 2087619, 1939205, -2153326, -2160422, -2162801, 2042194, -1993459, 2120517,
                1504726, 1306844, -2050497, -2085995, -554515, -1976499, -1906465, 1836648, 2075222, -1885203, -2153631,
                -2101753, 2129668, -2036718
            }//9
    };

    std::ofstream readfasta("read.fasta");
    HaplotypeConsensus haplotypeConsensus(ws, emptyLDG, emptyLDG, lines[5]);
    haplotypeConsensus.use_long_reads_from_file("reads_in_iline_5.fasta");

    std::string consensus;

    auto read = haplotypeConsensus.read_seqs.begin();

    for (int i = 0; read != haplotypeConsensus.read_seqs.end(); i++, ++read) {
        sglib::OutputLog() << "Processing read " << read->first << std::endl;
        readfasta << ">" << read->first << std::endl;
        readfasta << read->second << std::endl;
        ws.long_read_mappers[0].read_paths[read->first] = ws.long_read_mappers[0].create_read_path(read->first, true,
                                                                                                   read->second);
        haplotypeConsensus.oriented_read_paths.resize(
                std::max(read->first+1, (uint64_t) haplotypeConsensus.oriented_read_paths.size()));
        haplotypeConsensus.orient_read_path(read->first);
    }

    sglib::OutputLog() << "Done processing reads!" << std::endl << std::endl;
    sglib::OutputLog() << "Oriented read paths size: " << haplotypeConsensus.oriented_read_paths.size() << std::endl;

    haplotypeConsensus.build_line_path();
    consensus = haplotypeConsensus.consensus_sequence();
    std::ofstream consensusfasta("consensus5.fasta");
    consensusfasta << ">consensus" << std::endl;
    consensusfasta << consensus << std::endl;

    return 0;
}