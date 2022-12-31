import 'package:flutter/material.dart';
import 'package:flutter_blue_plus/flutter_blue_plus.dart';
import 'package:url_launcher/url_launcher.dart';

void main() {
  runApp(MyApp()); // 메인 페이지를 runApp에 넣어주면 됨
}
// StatelessWidget은 변화지 않는 화면을 작업할 때 사용.
// 변화는 화면을 작업 하고싶을 경우에는 StatefulWidget을 사용.
class MyApp extends StatelessWidget {
  const MyApp({Key? key}) : super(key:key);

  @override
  Widget build(BuildContext context) {
    return const MaterialApp(
        home: NavigationExample()
    );
    // return MaterialApp() -> Material 디자인 테마를 사용
  }
}
class NavigationExample extends StatefulWidget {
  const NavigationExample({Key? key}) : super(key: key);

  @override
  State<NavigationExample> createState() => _NavigationExampleState();
}



class _NavigationExampleState extends State<NavigationExample> {
  int currentPageIndex = 0;
  FlutterBluePlus flutterBlue = FlutterBluePlus.instance;
  final Uri _url = Uri.parse('https://portal.jnu.ac.kr');
  bool isSwitched = false;

  var backcolor = Colors.black;

  void _launchUrl() async {
    if (!await launchUrl(_url)) throw 'Could not launch $_url';
  }

  @override
  void initState() {
    super.initState();
    flutterBlue = FlutterBluePlus.instance;
  }

  void _startScan() {
    setState(() {});
    flutterBlue.startScan(timeout: Duration(seconds: 10));
    var subscription = flutterBlue.scanResults.listen((results) {
      // do something with scan results
      for (ScanResult r in results) {
        print('Device Name : ${r.device.name} // Device ID : ${r.device.id} // Device rssi: ${r.rssi}');
      }
    });
  }

  @override
  Widget build(BuildContext context) {
    return Scaffold(
      appBar: AppBar(
        centerTitle: true,
        title: Text('Jayou Haggi'),
        backgroundColor: Colors.lightGreen,
        actions: <Widget>[
          IconButton(
              onPressed: (){
                ScaffoldMessenger.of(context).showSnackBar(
                    const SnackBar(content: Text('자유학기 컨텐츠')));
              },
              icon: const Icon(Icons.abc),
              tooltip: 'show message'
          ),
        ],
      ),
      bottomNavigationBar: NavigationBar(
        onDestinationSelected: (int index) {
          setState(() {
            currentPageIndex = index;
          });
        },
        selectedIndex: currentPageIndex,
        destinations: const <Widget>[
          NavigationDestination(icon: Icon(Icons.home), label: 'home',),
          NavigationDestination(icon: Icon(Icons.accessibility), label: '개인 페이지',),
          NavigationDestination(icon: Icon(Icons.bluetooth), label: '블루투스 기기 검색',),
        ],
      ),
      body: <Widget>[
        Column(
          mainAxisAlignment: MainAxisAlignment.spaceEvenly,
          crossAxisAlignment: CrossAxisAlignment.stretch,
          children:[
            Container(
              height: 150,
              margin: EdgeInsets.all(10.0),
              child: Row(
                mainAxisAlignment: MainAxisAlignment.center,
                children: [
                  Container(
                    width: 360,
                    child: TextButton(
                      child: Image.asset('assets/jnu.jpg'),
                      onPressed: _launchUrl,
                    ),
                  ),
                ],
              ),
            ),
            Container(
              child: Row(
                mainAxisAlignment: MainAxisAlignment.spaceEvenly,
                children: <Widget> [
                  Container(
                    width: 300,
                    height: 150,
                    color: backcolor,
                  ),
                  Container(
                    margin: EdgeInsets.fromLTRB(0, 0, 10, 0),
                    child: Switch(
                        value: isSwitched,
                        onChanged: (value) {
                          setState((){
                            isSwitched = value;
                            if (isSwitched == true) {
                              backcolor = Colors.blueAccent;
                            } else {
                              backcolor = Colors.black;
                            }
                          });
                        }
                      ),
                  ),
                ],
              ),
            ),
            Container(
              height: 150,
              margin: EdgeInsets.all(10.0),
              child: ElevatedButton(
                style: ElevatedButton.styleFrom(primary: Colors.lightBlue),
                onPressed: () {
                  Navigator.push(
                    context,
                    MaterialPageRoute(builder: (context) => const SecondRoute())
                  );
                },
                child: Text('Function 실험용 버튼'),
              ),
            ),
          ],
        ),
        /*장치검색, 검색중지*/
        Container(
            alignment: Alignment.bottomCenter,
            color: Colors.white,
            child: Column(
              mainAxisAlignment: MainAxisAlignment.start,
              children: [
                Row(
                  children: <Widget> [
                    Container(
                      margin: EdgeInsets.fromLTRB(10, 5, 5, 0),
                      child: ClipRRect(
                        borderRadius: BorderRadius.circular(10),
                        child: Image.asset(
                          'assets/KakaoTalk_20220610_170941809.jpg',
                          width: 150,
                        ),
                      ),
                    ),
                    Expanded(
                      child: Container(
                        margin: EdgeInsets.fromLTRB(10, 10, 10, 10),
                        child: Column(
                          crossAxisAlignment: CrossAxisAlignment.start,
                          children: [
                            Text('이름: 김한진'),
                            Text("나이: 25"),
                            Text("키: 177"),
                            Text("몸무게: 71"),
                            Text("주소: 광주광역시"),
                            Text("기타 등등"),
                              ]
                        ),
                      ),
                    )
                  ],
                ),
                Expanded(
                    child: Column(
                      mainAxisAlignment: MainAxisAlignment.start,
                      children: [
                        Container(
                          margin: EdgeInsets.fromLTRB(10,7,10,0),
                          decoration: BoxDecoration(
                            borderRadius: BorderRadius.circular(10),
                            color: Colors.grey
                          ),
                          height: 100,
                          child: Row(
                            mainAxisAlignment: MainAxisAlignment.spaceEvenly,
                            children: [
                              Column(
                                mainAxisAlignment: MainAxisAlignment.center,
                                children: [
                                  Text('일일 활동'),
                                  Text('   '),
                                  Row(
                                    children: [
                                      Icon(Icons.account_tree),
                                      Text('      '),
                                      Icon(Icons.ac_unit_outlined),
                                      Text('      '),
                                      Icon(Icons.accessibility_sharp)
                                    ],
                                  ),
                                  Row(
                                    children: [
                                      Text('1200'),
                                      Text('      '),
                                      Text('20'),
                                      Text('      '),
                                      Text('100')
                                    ],
                                  )
                                ],
                              ),
                              Container(
                                width: 50,
                                child: Icon(Icons.access_time),
                              )
                            ],
                          ),
                        ),
                        Container(
                          margin: EdgeInsets.fromLTRB(10,7,10,0),
                          decoration: BoxDecoration(
                              borderRadius: BorderRadius.circular(10),
                              color: Colors.grey
                          ),
                          height: 100,
                          child: Row(
                            mainAxisAlignment: MainAxisAlignment.spaceEvenly,
                            children: [
                              Column(
                                mainAxisAlignment: MainAxisAlignment.center,
                                children: [
                                  Text('음식'),
                                  Text('   '),
                                  Row(
                                    children: [
                                      Text('200'),
                                      Text('     '),
                                      Text('/  1900kcal')
                                    ],
                                  ),
                                ],
                              ),
                              Container(
                                width: 60,
                                child: ElevatedButton(
                                  onPressed: () {},
                                  child: Text('입력'),),
                              )
                            ],
                          ),
                        ),
                        Container(
                          margin: EdgeInsets.fromLTRB(10,7,10,0),
                          decoration: BoxDecoration(
                              borderRadius: BorderRadius.circular(10),
                              color: Colors.grey
                          ),
                          height: 100,
                          child: Row(
                            mainAxisAlignment: MainAxisAlignment.spaceEvenly,
                            children: [
                              Column(
                                mainAxisAlignment: MainAxisAlignment.center,
                                children: [
                                  Text('심박수'),
                                  Row(
                                    children: [
                                      Text('90', style: TextStyle(fontSize: 30),),
                                      Text('  bpm')
                                    ],
                                  ),
                                ],
                              ),
                            ],
                          ),
                        ),
                      ],
                    ),
                ),
              ],
            ),
        ),
        Scaffold(
          body: Center(
            child: Column(
              mainAxisAlignment: MainAxisAlignment.center,
              children: <Widget> [],
            ),
          ),
          floatingActionButton: FloatingActionButton(
            onPressed: _startScan,
            tooltip: 'increment',
            child: Icon(Icons.bluetooth),
          ),
        ),
      ][currentPageIndex],
    );
  }
}
class SecondRoute extends StatelessWidget {
  const SecondRoute({super.key});

  @override
  Widget build(BuildContext context) {
    return Scaffold(
      body: Center(
        child: ElevatedButton(
          onPressed: () {
            Navigator.pop(context);
          },
          child: const Text("뒤로 가기"),
        ),
      ),
    );
  }
}