#!/usr/bin/python3
from random import randint
from math import log, sqrt
import time

# CPU low, MEM intensive

MAX_TRIES = 15


def sieve(n):
    try:
        psieve = [True] * n
        psieve[0] = psieve[1] = False  # 0 and 1 are not prime!
    except (MemoryError, OverflowError) as e:
        # raise Exception('Not enough memory to create the sieve!')
        n = n // 2
        print('Not enough memory to create the sieve! Reducing by half... Testing now with %d.' % n)
        return sieve(n)
    except IndexError as e:
        raise Exception('The number is too small to create the sieve!')

    primes = list()
    for (i, p) in enumerate(psieve):
        if p:  # if p is marked as prime
            primes.append(i)

            for j in range(i * i, n, i):
                psieve[j] = False  # remove all composites
    return primes


def gcd(m, n):
    while n != 0:
        m, n = n, m % n
    return m


def smooth_rho(n, primes):
    a = randint(2, n - 1)
    d = gcd(a, n)

    if d >= 2:
        return d

    for q in primes:
        l = log(n) // log(q)
        t = int(pow(q, l))
        a = pow(a, t, n)

    d = gcd(a - 1, n)

    return None if d == 1 or d == n else d


def find_factors(n, bound, tries=0):
    if n % 2 == 0:
        return 2, n // 2

    if tries > MAX_TRIES or bound > int(1 + sqrt(n)):
        #return 'No solution found using this bound.'
        return -1, -1

    primes = sieve(bound)

    r = smooth_rho(n, primes)
    f = n if r == 1 else r

    if r:
        return f, n // r
    else:
        bound *= 2
        #print('Increasing bound to', bound)
        return find_factors(n, bound, tries + 1)


if __name__ == '__main__':
    # p = randint(1, 2 ** 512)

    primes128 = [175120399228547846888456991422137905271, 35244174869138595766891618528472617627,
                 48333053950913730667611297492172966883, 267090296098723649371372079022126248161,
                 118453697632928227185360224459550204077, 195193436617336812952419535015699137637,
                 250937939158521665039181319524520106101, 278119091008525468679249984037864673841,
                 34891089573282708143478547396064064223, 197323530232852371312369746732591064827,
                 302078778715324364079338634276794360411, 66119965553181872836704158575992234313,
                 129903898072054084144016744719765296739, 140183990455859180343746800787622683233,
                 124690570437888380393695933578195039511, 20922036031177005386821330303691802343,
                 265749097497625777668593909381436791549, 1585877954466250064009619009913436203,
                 69223687829615033586804796110040130819, 201420733574071303955086719689417120413,
                 119467350619837503062006769459692977379, 243278481509601389448673216677201605299,
                 175623244331224589704743298745577914011, 228295726348180586272775652459490632635,
                 154389965900948591939526769839672790537, 183098177406452502925312700480373207547,
                 168497467037655466444508299295283974669, 307569403922904512865758297833966688917,
                 116692915686744231054806696697991688111, 136308016456696544131280852221353506153,
                 144464304873283791578400222799300004175, 28107318757437074658653670749738489405,
                 169426918302190818906859004094446621765, 193464348228454986762011307083414325185,
                 82353225050568614030107655572293222265, 246757809118530470110537523564892540395,
                 13489250576212515363541212504997327735, 270115781247430004585595378903218740595,
                 317958974102048441534356321576544708815, 106201146991857175246924168295138422915,
                 76087895074125393891670533036688318885, 330943842713466562465769765936218619755,
                 84615456425746900050890102923321848355, 185419019741466591923920010051165433865,
                 126773721129410087583167506895717604545, 9952641331915860561893816618877567755,
                 202842722466994589710540425697377286825, 278679376721486414894765804951069145365,
                 14802997725670076077805445098680091285, 233947908074071624097283218823329614675,
                 69974997296904050351639615080221493581, 310372879220072195786521381717950335037,
                 184635961444398376324192180818010556173, 170561351537293936584221176863915481403,
                 335091014454072490734451632488724095067, 270106047993349327121243496304986014181,
                 325121590089601414949397847766775304765, 235416726595463312634380349030191773775,
                 116602791033654997996699883429770713431, 88799405875253582140648849066356896905,
                 276685663279221964243156517129504423047, 260767117597044056595362259157254980323,
                 282257304862817053338565588366540633463, 332808284982460434895906321415598401429,
                 208093250377509198223838167467494150251, 275137695073878416664995611148057882439,
                 289580092196974914621689933312234044563, 16904698929602374627807749589450435665,
                 70930534781132413856038156896120255837, 135518828708391521579563131163435637685,
                 320643685134635070845969629689801386979, 76631144828097964098872353110661885665,
                 89463054422417761211843542642794881071, 133693059879703258553726187778664228169,
                 237812536467925333232844540523017751989, 335301107331811126431099902812435563611,
                 70902140475649942935511442705661431411, 147247646520807144078616823522859914149,
                 183000416444750885926064015719292170771, 183184338661231290295118985790535235525,
                 193667960621977349988521733081288602625, 222745343753038762660732071947563559919,
                 299033823683534210280253907129263499339, 485670170571956031721248064910607981,
                 3432164708076276110431832866460062599, 59256789692980160691624886475591445861,
                 63207415755041259253585795838940284049, 211909660886938714565712701863678062423,
                 220514623736472985103189384867711296311, 295559940224188918072594154467584810611,
                 308610470196113520632210715959039301239, 319700538072381326112161080537297855233,
                 32029559317965818601438484755201768813, 62494660749740024241025053169810004467,
                 64188525503393795897564176439371266573, 101437184301201157491420529456084197511,
                 273909789142190371260807445768865305291, 277407752498958233532929669804432921427,
                 147432365803156432856752081406983463241, 268200004857228315540149541008467670925]
    primes256 = [104719802762381912039208921868056844202718606698699371733963822856819421644967,
                 87284775599922252854173389331031823480309327938706104528688277212807448818173,
                 105359025658534907028220924351401714625261277413486767163373122866009445998593,
                 27512161426807114034151418592568159227183493100025344857129003297910113382799,
                 12432238490034514303265762210242373345321110965388236763711044337672967381737,
                 108522889293525336938684162998067024394018600701671805038455886663409539230263,
                 114480297316797915134262700236215164892791622152837215174432070427333400612481,
                 14818954088536007508519299979142982879679303666460180532218523273968473393197,
                 19962680072948101125825390522331074241461474044318698775705699243882742412057,
                 108746808525351060713022249092381272290976002096248954252526851005132206113971,
                 19467676898496091311044639903815498018587912259296715577427047660401410898873,
                 113325889372251743694489318928610285544633601230157863723807181339269679550171,
                 84715647122109906310826397841320371857928183409258502518089994073128694259617,
                 75096497439372458662193833060918583997611415480563231549209466046393345783113,
                 24800130526692521783273217524877579746709617446266957833608463857086933936359,
                 72395480271935102608498549571797057442537844770525733205601703124416039804151,
                 50275323910183797294050202110303055731076089533872096047815545139660947190657,
                 86893527126955781031220957861134716456953643948597990493900725833975863291379,
                 44427929182795530249399407506378792076640362772785741679041302191981411992777,
                 46300147893326657529660416883438637431131079032779575044529320677758896986103,
                 92239203038686861935846331939486221113474270817087800511915335320778587080847,
                 31470780473038591730655014950987598428553463339919352067078732048306916111607,
                 4435102415851265209952613896217197565479460450772528547863114613933821106327,
                 89091593665949220312324665712383116416986221856956278265146242314324650003943,
                 100047409087076394859933453556431209952214161767114552407020648438605605226833,
                 102878349753083477401932880024218690491620730847763134320033173076238921100973,
                 16787506313130862499226109639119080651306228056739263769826369871543608247709,
                 97676111825525356903846651698208421838970180612064960766708085766062193806123,
                 31700119697367123280008421933233753514303021477194513903377396136031130053963,
                 25008287202670808323567730971619412133059659217419071210688649207001360458883,
                 61236732854596897246160032145778059213703015241920688700963654636861107074983,
                 24959741237633420232762882511572386651207026561062179939884838834500083538169,
                 18991213446916390408410060093179066389860343014192733450910847830476341685769,
                 57759705118767493482621224128891749349329511557044189976723634402790884734523,
                 108513710649820811710498003165304484252262561772984785306921614123300981546931,
                 16736096961222077976767498555773732943286323739838088101666002569696467985903,
                 81218925777439164950981156713362927216977394514951054519234758847586865521789,
                 101147547872721129674753486585140186684665065435220232229381422638500591798335,
                 55436042209959939389547749148267209034505663972376599507226182839495019755575,
                 20709307321717105814311233997708316781901468051125763804854487332147075754165,
                 114269984751594828289304619535830087381871794805440064379591110494723361693235,
                 45635463143746123890070154148000154797918902251900510182950267681281309318825,
                 26574345733024469635595724811631175894439335888523825231538659800793438254145,
                 98707142070957779729502932044070398580197741440809870664727196959713121933835,
                 62360413641876958380568792351822183728094243735061968562311763890056586454935,
                 77828424936270456006659664954681543092406608345703697607212691073272942378755,
                 83261983938159634305360950637454053058799374737035477068477321752481396312121,
                 68327806421898589341189154360269977144118000032162215548396218656697127147305,
                 11157695570649416715220898826122721017336640914203862713122681246993837275077,
                 97011566502266201270181312347178952095599514968587120739746511648598077682041,
                 105093097720703975358749445817524269106839080945308786293820640838103298233065,
                 27951938608906250768721952326638093443592565655126561768781633311381344377945,
                 64622496188442161057182195522852959660434681919364264917142915494074462439433,
                 13203696853326809082063425460190208005957618139091699063747386881449527698895,
                 30643725465148189523381163411638231723872854772368713719220004667587991816913,
                 53735288970390104763139004818895068199849773424742231325216972802174605357931,
                 57919854749579060379453169442312740280974397122499846034163366833269201636609,
                 19878967925299363913691397618049869324695823353776458625487297763892615944105,
                 100981945499678949163130309482855180506963092188587891605687097413086116000867,
                 74869937407432760072504886492280509428501011970499847660415944581149200160835,
                 16602870209027105204703330599278592911440930296575992896464987920311360506527,
                 95464622069459096171528824133523461014416524622882426868786189327623434668227,
                 27064069831021020466544958881309001310205527637395645343800119527111235740923,
                 42827666082012872025559168911369623232493864667827808399018287726356956536643,
                 55048144715158215939466453554318982670410034986630543617371022391238489949099,
                 63813771019720559323572573475948356766584350833084529007580719985908215798829,
                 69153569008007802503471079004435483564353669860644177123777648395787228972041,
                 32398545877780884462462701064743481028816121082402710352350676446802587176455,
                 1925403017322192669118485022690898477252275775251847921555046950984515062745,
                 3389442495097181201842597592144898055039304233970732829791359492568372941341,
                 6902610959649186812567120781721680874780354062253363161690859139907525836883,
                 11245053027956132101125321800317955960409193639624649943029859309393435680475,
                 13455689126174303930628275492733298152282028688266343887022282048368126409871,
                 30649436144676824318456761729092415386129381655789445406073091262527321805633,
                 33434548252489473080919881770971449952026122221110184096147669969252689846349,
                 35562310688991314941085225489322509053103306388099679030374922719137006131143,
                 38038203661336113257445207193691928434558639787037769000790270630935687113629,
                 78098063587365910982256349083608903026293004963509667569476519916148532252103,
                 90780116097150938671161691778971561458843031786835889430343632765770993297243,
                 95803791234833099969122091011064488714673766123059708190348447424598699570343,
                 105661963751904584073868802590204469221612765375656086455856419523591681088357,
                 108989907511970108572917948653649670490709592863954706923650663525672347071227,
                 112676275318238626722781327150388838642073629426366341321157355921309628695533,
                 16542882936380246952618429325548914089747295212896773179608378824050552431167,
                 50906072632548030583810820599462066816361191135869625245352782152359739920881,
                 57719791604083843972379524809323056733716970007499082594239511864002577441517,
                 67933280214237296542215049660150910517869535930150651052501679250972386933091,
                 72973343624894346317151228235675724749769332092552075444635695165649280747949,
                 77557113050614327229384432076336533286421568631022477582391542831955457986133,
                 93300371954938480124412237272833608967872698925783153475709702043568653244635,
                 101228219946727387224239298785762275013719216945208380553014448472957944791713,
                 107410958565959202744390531474634638684079781818364119408298214215184028343523,
                 114010526990660219947796123520266408662815205140541809571378492429240718081281,
                 13046098600387942205237366796948221748934879203700800470042985619679099807133,
                 22624166035122230078729551042935246959872003635816098899455893352776635388961,
                 42503680158145375576333131105069778451023704189808463469327591560211314764329,
                 72130696849225938196910438889980376404231078214442037893045587608086077001503,
                 75344259581306715324370347765821912437734782447456629872322893204624889325859,
                 106996548023607834073900582087233730663908554890584512900206156280316515248383,
                 113308701910119931509489914452025447668454635109148558508642780071520865423877]

    for p in primes256:
        smooth_bound = int(log(p, 2) + 1) if p % 2 == 1 else 2
        #print('prime: %d bound: %d' % (p, smooth_bound))
        #print(find_factors(p, smooth_bound))
        ini = time.time()
        f1, f2 = find_factors(p, smooth_bound)
        elapsed_time = time.time() - ini

        print('%d ; %d, %d ; %4.8f' % (p, f1, f2, elapsed_time))