A Nemzeti Szívinfarktusregiszter (NSZR) interaktív lekérdező felülete
================
Ferenci Tamás (<https://www.medstat.hu/>)

A Nemzeti Szívinfarktusregiszter (NSZR) interaktív lekérdező felülete a
<https://nszr.gokvi.hu/lekerdezo/> címen érhető el.

Jelen oldal célja kettős: egyrészt bemutatni a lekérdező felület
háttérét, a létrehozása mögötti motivációt, másrészt rövid magyarázattal
szolgálni az elérhető funkciókhoz, elemzésekhez (kitérve azok korlátaira
is). E mellett itt érhető el – a transzparencia és a nyílt tudomány
jegyében – a lekérdező oldal hátterében lévő program teljes forráskódja
is, mely megkönnyíti, hogy az esetleges hibái kiderüljenek, illetve
elősegíti, hogy jobb megoldások, továbbfejlesztési lehetőségek
merüljenek fel.

## Az NSZR lekérdező felületének létrehozása mögötti motiváció

Mind a tudományos megismerés, mind a hatékony népegészségügyi munka
szempontjából alapvető fontosságú, hogy jól értsük a betegségek
epidemiológiáját, tehát tudjunk válaszolni arra, hogy kik, mikor, hol,
és miért betegek. (Sok esetben ez utóbbi, talán legfontosabbnak tűnő
kérdés megválaszolásában pont az a nagyon fontos információ, ha az
előbbiekre van jó válaszunk.)

Az egészségügyi ellátás jellemzése, az ellátórendszer javítási
lehetőségek feltárása szempontjából két kérdés nagyon fontos: egyrészt,
hogy jó helyzetképünk legyen az ellátás aktuális jellemzőiről (ideértve
összefüggéseit, például a betegjellemzőkkel, időponttal, helyszínnel
stb.), másrészt, hogy legyen információnk a kimenetről, például a beteg
túléléséről. Ez utóbbi azért fontos, mert egyszerre mutatója a beteg
állapotának és az ellátás minőségének, de ha statisztikai eszközökkel
szét tudjuk választani a kettőt, akkor kulcsfontosságú információkat
nyerhetünk az ellátás jóságáról, ami – különösen, ha összefüggésbe
hozzuk egyéb jellemzőkkel, például helyével vagy idejével – nagyban
segítheti a javítási lehetőségek feltárását is.

Mindezek megválaszolásához azonban információra van szükségünk. Sokféle,
lehetőség szerint jó minőségű információra. Hogyan tudjuk ezeket
begyűjteni?

Nagyon általános szinten tekintve három lehetőség jöhet szóba:

- Az egyik megoldási lehetőség ad hoc vizsgálatok szervezése. (Például
  egy mintavétellel történő felmérés – jó esetben véletlenszerűen a
  populációból, rosszabb esetben ún. kényelmi mintaként, például
  egyetlen kórház adatainak feldolgozásával.) Ez kevesebb erőforrást
  igényel, de csak egy pillanatfelvételt ad, és semmiképp nem
  teljeskörű, kényelmi minta esetén pedig erősen kérdéses is az
  általánosíthatósága.
- A másik lehetőség az ún. adminisztratív/finanszírozási adatok
  felhasználása. Az alapötlet, hogy a kórházak amúgy is jelentenek
  finanszírozási célból adatokat – miért ne lehetne ezt felhasználni
  epidemiológiai célokra is? Csakugyan, ha valaki egészségügyi
  ellátásban részesül, akkor keletkezik róla egy adatsor, benne a
  nemével, életkorával, lakhelyével, betegségével, az elvégzett
  beavatkozással; ebből tényleg kiolvasható lehet egy sor hasznos
  információ. Ez nagyon csábítóan hangzik, hiszen az erőforrás-igénye
  csekély (ezeket az adatokat amúgy is begyűjtik), de mégis teljeskörű
  és folyamatosan frissülő az adatbázis, legalábbis a közfinanszírozott
  ellátásokra vonatkozóan. Bár ez eddig jól hangzik, de a módszernek
  vannak hátrányai is: egyrészt az adatminőség (ezeket a jelentéseket a
  kórházak rutinszerűen „optimalizálják” finanszírozási szempontból),
  másrészt a klinikai adatok hiánya (azt tudjuk, hogy valakit
  megröntgeneztek, de azt nem tudjuk, hogy mi volt a röntgenképen; azt
  tudjuk, hogy az alany hány éves, de azt nem tudjuk, hogy dohányzik-e
  stb.). Ezzel együtt is, ma már egyre több ilyen vizsgálat készül; egy
  példa tisztán akadémiai célokat szolgáló ilyen kutatásra a
  [HUNVASCDATA](https://hunvascdata.hu/)-projekt.
- Végezetül a harmadik lehetőség a betegségregiszterek használata. A
  regiszter definíció szerint azt jelenti, hogy valamely megbetegedés
  előfordulásáról a teljesség igényével készül gyűjtés, klinikai
  adatokkal együtt. (Tipikusan jogszabály írja el a kötelező jelentést
  az egészségügyi ellátóknak.) Ez látszólag az ideális megoldás:
  teljeskörű, folyamatos, validálható adatminőségű, részletgazdag
  klinikai adatokat is tartalmazhat, egyetlen problémája van: az, hogy
  hatalmas az erőforrásigénye. Nem csak „forintban” értve, hanem az
  adatszolgáltatói teherre nézve is, hiszen ez azt is jelenti, hogy az
  észlelő orvosoknak minden beteg után egy plusz jelentést is ki kell
  tölteniük, és feltölteni a regiszterbe.

A fentiekből is érzékelhető, hogy a regiszter használata, bár sok
szempontból a legjobb lehetőség, jellemzően csak a népegészségügyileg
legfontosabb betegségeknél jöhet szóba.

A szívinfarktus azonban egyértelműen ilyen, ezért – néhány korábbi,
korlátozottabb kísérlet után – 2010-ben elindult mai formájában a
Nemzeti Szívinfarktusregiszter (NSZR), amelybe 2014. január 1. óta
törvényileg kötelező jelenteni minden észlelt szívinfarktusos esetet,
részletes klinikai adatokkal együtt. A regiszer azóta is töretlenül
működik, végezve a teljeskörűség igényével az adatgyűjtést (és a
gyűjtött adatok validációját is).

Az adatok gyűjtése azonban csak a feladat egyik része: fel is kell
azokat használni. Az NSZR rendszeres jelentéseket készít az irányító
egészségügyi hatóságoknak, valamint, ami a mi szempontukból még
fontosabb, ezekből az adatokból számos tudományos publikáció is készült
az évek során, ami ezt a hasznosítási célt részben megvalósítja.
Részben, de nem tökéletesen: az egyik probléma, hogy a tudományos
publikációk mindig csak egy pillanatképet adnak, a tartalmuk nem
frissül, a másik, hogy mindig csak egy adott témáról szólnak, amit
ráadásul a kutatók határoztak meg, az olvasó nem tudja megváltoztatni.
Szerencsére a mai informatikai lehetőségek mellett mindkét hiányosság
pótolható: létrehozható egy olyan webes, interaktív felület, mely
folyamatosan frissülő adatokat tartalmaz, és melyen az érdeklődő saját
maga választhatja meg az elemzés tárgyát, beállításait. Fontos
hangsúlyozni, hogy ez nem váltja ki a tudományos publikációkat,
ellenkezőleg, a két dolog szerepe egymást kiegészítő jellegű: a
publikációk ugyanis a fenti nehézségeik mellett egyedi adatokon
alapulnak, amik nyilvánosan nem közölhetőek, de finomabb elemzést
tesznek lehetőve, míg a lekérdező felület a fenti előnyök mellett
mutathatja be a kevésbé részletezett, aggregált adatokat.

Egy ilyen interaktív lekérdező felület létrehozása és bárki számára
hozzáférhető működtetése számos előnnyel és népegészségügyi haszonnal
bír, melyek révén szolgálja az ország egészségét.

Az első a közvélemény informálása: a lakosság teljes joggal várja el,
hogy az ő pénzéből az ő érdekében gyűjtött adatok megismerhesse. Mindez
azonban nem öncél: az ilyen adatok közlése ismeretterjesztés is, mely
elősegíti, hogy a lakosság jobban megértsen egy betegséget, elmélyíti az
egészségtudatosságot, kapcsolatot épít betegek és az egészségügyi
ellátórendszer között. Fontos, hogy a weboldal nem csak eredményeket
közöl, hanem minden ponthoz magyarázot is ad, ami igyekszik az
elemzéseket minél több érdelődő számára hozzáférhetővé tenni, az
eredményeket kontextusba helyezni, valamint megerősíti az
ismeretterjesztő jelleget is. Itt érdemes hozzátenni, hogy az ilyen
adatközlés csapdákat is rejthet magában, hiszen az adatok több esetben
félreérthetőek lehetnek, ha kellő háttérismeret nélkül tekintenek rájuk,
azonban egy jó felület ezt a problémát pont hogy csökkenti, mert
magyarázatok révén megadhatja ezeket a háttérismereteket (amelyek révén
ráadásul saját magán is túlmutató felvilágosítást is végezhet).
Valójában ráadásul többről is szó van, mint a lakosság informálásáról.
Az adatokat ugyanúgy használhatják kutatók, elemzők, orvosok,
egészségügyi dolgozók, így egyrészt maga is lehet kutatások adatforrása,
ideértve akár a tájékozódást is például az aktuális népegészségügyi
helyzetről, másrészt megtermékenyíthetőleg hathat más kutatások
elvégzésére is.

A második szempont a bizalom megteremtése: a népegészségügyi munka egyik
alappillére a bizalom megteremtése és fenntartása, ennek pedig az egyik
legjobb eszköze a nyílt és őszinte kommunikáció, a transzparencia. Az
NSZR lekérdező felülete erre kettős példát ad: nem csak az adatokat
teszi elérhetővé, de a lekérdezést megvalósító felület, az annak
hátterében lévő számításokat végző programkódokat is nyilvánosságra
hozza. (Az adatok esetében természetesen az adatvédelmi szabályoknak
megfelelés érdekében kizárólag aggregált adatokat közölve.)

A lekérdező felület létrehozatalának harmadik célja a példaadás:
jogszabály szerint Magyarországon több mint egy tucat regiszter kell,
hogy [működjön](https://njt.hu/jogszabaly/2018-49-20-5H), azonban ezek
legnagyobb részének nincs értelmezhető, ténylegesen teljeskörű,
folyamatosan frissülő, pláne kívülről is látható, bárki számára elérhető
aktivitása (van olyan regiszter, aminek a nevére rákeresve interneten az
egyetlen találat a törvény szövege…). Elfogadva a fenti érveket, ez egy
hiányosság – azonban az NSZR példája mutatja, hogy egy javítható
hiányosság! Az NSZR interaktív lekérdező felülete példát ad arra, hogy
minimális erőforrásból is megvalósítható egy kulturált, bárki számára
elérhető, jól használható weboldal, mely a fenti célok szolgálata révén
segíti a betegeket, az érdeklődőeket, az egészségügyi dolgozókat, a
népegészségügyet – az ország lakosságának egészét.

## A lekérdező felületen elérhető elemzések

Az NSZR interaktív lekérdező felületén elérhető elemzések szorosan
tükrözik az előző pont elején felvázolt feladatokat, célokat. E
pillanatban négy különböző aspektusból vizsgálathatóak az NSZR adatai.

- Az Előfordulás pont teszi vizsgálhatóvá a talán legalapvetőbb kérdést:
  hány megbetegedés történik? Ez a ,,kik betegek?’’ kérdéskör
  kiindulópontja. Kiválasztható, hogy az összes szívinfarktust
  ábrázolja, az ST-elevációs infarktusokat (STEMI) vagy az
  ST-elevációval nem járó infarktusokat (NSTEMI). Az események abszolút
  száma (esetszáma) ugyan érdekes adat, de nem összehasonlítható sem
  időben, sem térben, hiszen eltér a különböző időpontok, illetve
  különböző megyék lélekszáma. Éppen ezért fontos a következő mutató, a
  nyers ráta, ami nem más, mint az esetszám osztva a megfelelő
  lélekszámmal. Ez azonban még mindig nem tökéletes, hiszen
  elképzelhető, hogy különböző időpontok, illetve különböző megyék
  között nem csak a lakosok száma, de összetétele is eltér! Ha például
  az egyik megyében a lakosok idősebbek, akkor egy ottani nagyobb nyers
  ráta nem feltétlenül jelent rosszabb helyzetet, hiszen lehet, hogy
  csak az idősebb – és emiatt nagyobb szívinfarktus-kockázatú – lakosok
  miatt van. Ez egy példa arra az általános jelenségre, amit magyarban
  is gyakran használt angol szóval
  [confounding-nak](https://tamas-ferenci.github.io/FerenciTamas_AzOrvosiMegismeresModszertanaEsAzOrvosiKutatasokKritikusErtekelese/)
  szoktak nevezni. A standardizált ráta használata ezt szűri ki; fontos
  hangsúlyozni, hogy csak az életkor vonatkozásában. (Elvileg az is
  felvethető, hogy ha valamelyik megyében vagy valamelyik időpontban
  több cukorbeteg van, akkor ennek a hatását is szűrjük ki, hiszen
  lehet, hogy csak emiatt nagyobb a ráta, és nem a helyzet rosszabb
  ténylegesen. Ezt azonban ha akarnánk sem tudnánk megtenni, hiszen arra
  vonatkozóan nincs megbízható adatunk, hogy hány cukorbeteg van az
  országban, pláne nem megye, hónap és életkor-csoport szerint
  lebontva.)
- A második pont, a Társbetegségek tovább finomítja a ,,kik betegek?’’
  kérdésre adott választ. E tekintetben ugyanis alapvető az a kérdés,
  hogy az infarktust elszenvedőknek milyen egyéb betegségeik,
  kórelőzményi adataik vannak (ideértve az egyéb – krónikus –
  betegségeket, amik fennálltak az infarktus pillanatában, valamint a
  kórelőzményben szereplő fontosabb akut eseményeket). Az ideális az
  lenne, ha meg tudnánk mondani, hogy e betegségek milyen kockázatot
  jelentenek, például, hogy a cukorbetegek körében mennyi az infarktus
  előfordulása (illetve, hogy ez mennyivel nagyobb, mint a nem
  cukorbetegek körében). Amint az előző pontban már volt róla szó, ezt
  sajnos nem tudjuk megtenni: azt tudjuk, hogy mennyi infarktusos volt
  cukorbeteg, de nem tudjuk, hogy mivel kellene leosztani, mert nem
  tudjuk, hogy hány cukorbeteg van, pláne nem megyénként, havonta,
  életkor-csoportonként. Ezért jobb híján a fenti információ
  fordítottját tudjuk megadni: nem azt, hogy a cukorbetegek körében
  milyen gyakori az infarktus, hanem azt, hogy az infarktusosok körében
  milyen gyakori a cukorbetegség. Ez érhető el ebben a pontban; összesen
  8 különböző társbetegségre, illetve kórelőzményre vonatkozóan.
- A harmadik pont, az Ellátás áttér az ellátási helyzetkép kérdésére.
  Jelenleg itt egyetlen ellátási jellemző vizsgálható: az, hogy milyen
  arányban történt a betegeknél katéteres érmegnyitás (PCI). Akár
  térbeli, akár időbeli összehasonlítást végzünk, felmerül pontosan
  ugyanaz a confounding-probléma, mint az előfordulás vizsgálatánál: a
  PCI elvégzésének valószínűsége kapcsolatban van vagy lehet a betegek
  életkorával, nemével, társbetegségeivel, márpedig az összehasonlított
  helyek és időpontok eltérhetnek e szerint. Éppen ezért az egyszerű
  arányon kívül mód van ún. korrigált arányt is választani mutatóként,
  mely statisztikai eszközökkel szűri ezek hatását (beállíthatóan csak a
  nem és az életkor, vagy nem, életkor és társbetegségek). Ez a mutató
  jobban összehasonlítható lehet időben, illetve térben.
- Végezetül, a Túlélés pont lehetővé teszi az egyik legfontosabb
  kimeneti jellemző, a túlélés vizsgálatát. E pillanatban a 30 napos és
  az 1 éves túlélés elemezhető; mindkét esetben választható a nyers
  arányon túl az előzőekben látottaknak megfelelő korrigált arány is.
  Itt is elmondható, hogy – a confounding szűrése révén – a korrigált
  arányok jobban összehasonlíthatóak, akár térben, akár időben
  vizsgálódunk.

Az összes fenti esetben igaz, hogy a választott kérdés vizsgálható
időben (beállíthatóan éves vagy havi szinten), elemezve a „mikor
betegek?” kérdést, valamint térben (megyei szinten), elemezve „hol
betegek?” kérdést. Az időbeli adatok lebonthatóak különböző jellemzők
szerint, ami azt jelenti, hogy egyetlen görbe helyett többet kapunk – a
lebontás alapján szolgáló változó minden lehetséges értékéhez egyet.
(Például nemi lebontás választása esetén egyetlen görbe helyett két
görbénk lesz, egy férfi és egy női.) A térbeli adatok, tehát a térképes
ábrázolás esetén pedig szűkíthetünk, ami azt jelenti, hogy
kiválaszthatunk egy vagy több változót, és azoknak megadhatunk egy
értéket – ezt követően a térképen csak azok az eredmények fognak
szerepelni, ahol a vizsgált változó a beállított értékű. (Például csak a
férfiak adatai jelennek meg, vagy csak a cukorbetegeké; ezek egymással
kombinálhatóak is.)

Az eredmények minden esetben megjeleníthetőek ábraként és táblázatként,
mindkettő letölthető különböző formátumokban.

## A lekérdező felület forráskódja

Az NSZR lekérdező felületének háttérében egy [R
programnyelven](https://www.youtube.com/@FerenciTamas/playlists?view=50&sort=dd&shelf_id=2),
közelebbről [R Shiny](https://shiny.posit.co/) csomag segítségével
megírt alkalmazás áll. Az interaktív vizualizációkhoz a weboldal a
[highcharts](https://www.highcharts.com/) könyvtárat használja, a
[higcharter](https://jkunst.com/highcharter/) csomagon keresztül. Az
adatkezelést a [data.table](https://rdatatable.gitlab.io/data.table/)
csomag segíti.

A transzparencia és a nyílt tudomány jegyében ennek teljes forráskódja
elérhető a következő linken:
[app.R](https://github.com/tamas-ferenci/NSZR-Lekerdezo/blob/main/app.R).
Ennek tartalma kettős: egyrészt benne van magának a weboldalnak a
forráskódja (ideértve a felületet, az eredmények előállítását és
megjelenítését), másrészt benne van a háttérben lévő statisztikai
számításoknak – különféle szűrések, lebontások, korrekciók – forráskódja
is.

A forráskód közzétételének a célja, hogy megkönnyítse, hogy az esetleges
hibái kiderüljenek, illetve elősegítse, hogy javítások,
továbbfejlesztések születhessenek hozzá.

------------------------------------------------------------------------

Írta: Ferenci Tamás. E leírás utoljára frissítve: 2024. október 2.
