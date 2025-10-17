7. Verifying that diagnosis

Look in one of your pi0_corrected_counts_<period>.json files:
    •   Do you see non-zero "contamination":{"+1":{"value":…}} entries?

--

No, the contamination is non zero, e.g.:

{
  "binning_meta": {"phi_bins": 12, "xB_bins": 8, "Q2_bins": 7, "t_bins": 6},
  "bins": {
    "(0,0,0,0)": {"N_data":{"helicity":{"+1":1,"-1":0},"total":1},"contamination":{"+1":{"value":0.00000000,"err":0.00000000},"-1":{"value":0.00000000,"err":0.00000000}},"N_corrected":{"helicity":{"+1":{"value":1.00000000,"err":1.00000000},"-1":{"value":0.00000000,"err":0.00000000}},"total":{"value":1.00000000,"err":1.00000000}}},
    "(0,0,0,1)": {"N_data":{"helicity":{"+1":3,"-1":1},"total":4},"contamination":{"+1":{"value":0.00000000,"err":0.00000000},"-1":{"value":0.00000000,"err":0.00000000}},"N_corrected":{"helicity":{"+1":{"value":3.00000000,"err":1.73205081},"-1":{"value":1.00000000,"err":1.00000000}},"total":{"value":4.00000000,"err":2.00000000}}},
    "(0,0,0,2)": {"N_data":{"helicity":{"+1":0,"-1":2},"total":2},"contamination":{"+1":{"value":0.00000000,"err":0.00000000},"-1":{"value":0.00000000,"err":0.00000000}},"N_corrected":{"helicity":{"+1":{"value":0.00000000,"err":0.00000000},"-1":{"value":2.00000000,"err":1.41421356}},"total":{"value":2.00000000,"err":1.41421356}}},
    "(0,0,0,3)": {"N_data":{"helicity":{"+1":2,"-1":1},"total":3},"contamination":{"+1":{"value":0.00000000,"err":0.00000000},"-1":{"value":0.00000000,"err":0.00000000}},"N_corrected":{"helicity":{"+1":{"value":2.00000000,"err":1.41421356},"-1":{"value":1.00000000,"err":1.00000000}},"total":{"value":3.00000000,"err":1.73205081}}},
    "(0,0,0,4)": {"N_data":{"helicity":{"+1":9,"-1":8},"total":17},"contamination":{"+1":{"value":0.05555556,"err":0.07579030},"-1":{"value":0.03125000,"err":0.04815948}},"N_corrected":{"helicity":{"+1":{"value":8.49999996,"err":2.91428472},"-1":{"value":7.75000000,"err":2.76699295}},"total":{"value":16.24999996,"err":4.01861984}}},
    "(0,0,0,5)": {"N_data":{"helicity":{"+1":10,"-1":8},"total":18},"contamination":{"+1":{"value":0.03571429,"err":0.05524594},"-1":{"value":0.10000000,"err":0.13964240}},"N_corrected":{"helicity":{"+1":{"value":9.64285710,"err":3.09898059},"-1":{"value":7.20000000,"err":2.77992806}},"total":{"value":16.84285710,"err":4.16313352}}},
    "(0,0,0,6)": {"N_data":{"helicity":{"+1":14,"-1":9},"total":23},"contamination":{"+1":{"value":0.00000000,"err":0.00000000},"-1":{"value":0.00000000,"err":0.00000000}},"N_corrected":{"helicity":{"+1":{"value":14.00000000,"err":3.74165739},"-1":{"value":9.00000000,"err":3.00000000}},"total":{"value":23.00000000,"err":4.79583152}}},
    "(0,0,0,7)": {"N_data":{"helicity":{"+1":3,"-1":4},"total":7},"contamination":{"+1":{"value":0.00000000,"err":0.00000000},"-1":{"value":0.00000000,"err":0.00000000}},"N_corrected":{"helicity":{"+1":{"value":3.00000000,"err":1.73205081},"-1":{"value":4.00000000,"err":2.00000000}},"total":{"value":7.00000000,"err":2.64575131}}},
    "(0,0,0,8)": {"N_data":{"helicity":{"+1":1,"-1":0},"total":1},"contamination":{"+1":{"value":0.00000000,"err":0.00000000},"-1":{"value":0.00000000,"err":0.00000000}},"N_corrected":{"helicity":{"+1":{"value":1.00000000,"err":1.00000000},"-1":{"value":0.00000000,"err":0.00000000}},"total":{"value":1.00000000,"err":1.00000000}}},
    "(0,0,0,10)": {"N_data":{"helicity":{"+1":3,"-1":3},"total":6},"contamination":{"+1":{"value":0.00000000,"err":0.00000000},"-1":{"value":0.00000000,"err":0.00000000}},"N_corrected":{"helicity":{"+1":{"value":3.00000000,"err":1.73205081},"-1":{"value":3.00000000,"err":1.73205081}},"total":{"value":6.00000000,"err":2.44948974}}},
    "(0,0,0,11)": {"N_data":{"helicity":{"+1":2,"-1":0},"total":2},"contamination":{"+1":{"value":0.00000000,"err":0.00000000},"-1":{"value":0.00000000,"err":0.00000000}},"N_corrected":{"helicity":{"+1":{"value":2.00000000,"err":1.41421356},"-1":{"value":0.00000000,"err":0.00000000}},"total":{"value":2.00000000,"err":1.41421356}}},
    "(0,0,1,0)": {"N_data":{"helicity":{"+1":7,"-1":4},"total":11},"contamination":{"+1":{"value":0.00000000,"err":0.00000000},"-1":{"value":0.00000000,"err":0.00000000}},"N_corrected":{"helicity":{"+1":{"value":7.00000000,"err":2.64575131},"-1":{"value":4.00000000,"err":2.00000000}},"total":{"value":11.00000000,"err":3.31662479}}},
    "(0,0,1,1)": {"N_data":{"helicity":{"+1":17,"-1":12},"total":29},"contamination":{"+1":{"value":0.08823529,"err":0.12137267},"-1":{"value":0.09090909,"err":0.13145433}},"N_corrected":{"helicity":{"+1":{"value":15.50000007,"err":4.28832205},"-1":{"value":10.90909092,"err":3.52217406}},"total":{"value":26.40909099,"err":5.54936177}}},
    "(0,0,1,2)": {"N_data":{"helicity":{"+1":15,"-1":7},"total":22},"contamination":{"+1":{"value":0.00000000,"err":0.00000000},"-1":{"value":0.00000000,"err":0.00000000}},"N_corrected":{"helicity":{"+1":{"value":15.00000000,"err":3.87298335},"-1":{"value":7.00000000,"err":2.64575131}},"total":{"value":22.00000000,"err":4.69041576}}},
    "(0,0,1,3)": {"N_data":{"helicity":{"+1":10,"-1":11},"total":21},"contamination":{"+1":{"value":0.00000000,"err":0.00000000},"-1":{"value":0.00000000,"err":0.00000000}},"N_corrected":{"helicity":{"+1":{"value":10.00000000,"err":3.16227766},"-1":{"value":11.00000000,"err":3.31662479}},"total":{"value":21.00000000,"err":4.58257569}}},
    "(0,0,1,4)": {"N_data":{"helicity":{"+1":56,"-1":47},"total":103},"contamination":{"+1":{"value":0.01846154,"err":0.01327171},"-1":{"value":0.01573427,"err":0.01178694}},"N_corrected":{"helicity":{"+1":{"value":54.96615376,"err":7.38266643},"-1":{"value":46.26048931,"err":6.77048853}},"total":{"value":101.22664307,"err":10.01714922}}},
    "(0,0,1,5)": {"N_data":{"helicity":{"+1":60,"-1":60},"total":120},"contamination":{"+1":{"value":0.08163265,"err":0.03904239},"-1":{"value":0.05952381,"err":0.02899620}},"N_corrected":{"helicity":{"+1":{"value":55.10204100,"err":7.48942087},"-1":{"value":56.42857140,"err":7.48976198}},"total":{"value":111.53061240,"err":10.59188177}}},
    "(0,0,1,6)": {"N_data":{"helicity":{"+1":40,"-1":38},"total":78},"contamination":{"+1":{"value":0.05664488,"err":0.03552731},"-1":{"value":0.03921569,"err":0.02564238}},"N_corrected":{"helicity":{"+1":{"value":37.73420480,"err":6.13320951},"-1":{"value":36.50980378,"err":6.00229309}},"total":{"value":74.24400858,"err":8.58159550}}},
    "(0,0,1,7)": {"N_data":{"helicity":{"+1":17,"-1":21},"total":38},"contamination":{"+1":{"value":0.19780220,"err":0.16781511},"-1":{"value":0.04285714,"err":0.04341474}},"N_corrected":{"helicity":{"+1":{"value":13.63736260,"err":4.36791193},"-1":{"value":20.10000006,"err":4.47993145}},"total":{"value":33.73736266,"err":6.25687146}}},
    "(0,0,1,8)": {"N_data":{"helicity":{"+1":4,"-1":4},"total":8},"contamination":{"+1":{"value":0.00000000,"err":0.00000000},"-1":{"value":0.00000000,"err":0.00000000}},"N_corrected":{"helicity":{"+1":{"value":4.00000000,"err":2.00000000},"-1":{"value":4.00000000,"err":2.00000000}},"total":{"value":8.00000000,"err":2.82842712}}},
    "(0,0,1,9)": {"N_data":{"helicity":{"+1":14,"-1":13},"total":27},"contamination":{"+1":{"value":0.08333333,"err":0.14632852},"-1":{"value":0.07692308,"err":0.13493200}},"N_corrected":{"helicity":{"+1":{"value":12.83333338,"err":3.99507797},"-1":{"value":11.99999996,"err":3.76215974}},"total":{"value":24.83333334,"err":5.48766744}}},
    "(0,0,1,10)": {"N_data":{"helicity":{"+1":15,"-1":14},"total":29},"contamination":{"+1":{"value":0.00000000,"err":0.00000000},"-1":{"value":0.00000000,"err":0.00000000}},"N_corrected":{"helicity":{"+1":{"value":15.00000000,"err":3.87298335},"-1":{"value":14.00000000,"err":3.74165739}},"total":{"value":29.00000000,"err":5.38516481}}},
    "(0,0,1,11)": {"N_data":{"helicity":{"+1":5,"-1":6},"total":11},"contamination":{"+1":{"value":0.00000000,"err":0.00000000},"-1":{"value":0.00000000,"err":0.00000000}},"N_corrected":{"helicity":{"+1":{"value":5.00000000,"err":2.23606798},"-1":{"value":6.00000000,"err":2.44948974}},"total":{"value":11.00000000,"err":3.31662479}}},
    "(0,0,2,0)": {"N_data":{"helicity":{"+1":32,"-1":38},"total":70},"contamination":{"+1":{"value":0.00000000,"err":0.00000000},"-1":{"value":0.00000000,"err":0.00000000}},"N_corrected":{"helicity":{"+1":{"value":32.00000000,"err":5.65685425},"-1":{"value":38.00000000,"err":6.16441400}},"total":{"value":70.00000000,"err":8.36660027}}},
    "(0,0,2,1)": {"N_data":{"helicity":{"+1":84,"-1":58},"total":142},"contamination":{"+1":{"value":0.00000000,"err":0.00000000},"-1":{"value":0.00000000,"err":0.00000000}},"N_corrected":{"helicity":{"+1":{"value":84.00000000,"err":9.16515139},"-1":{"value":58.00000000,"err":7.61577311}},"total":{"value":142.00000000,"err":11.91637529}}},
    "(0,0,2,2)": {"N_data":{"helicity":{"+1":57,"-1":46},"total":103},"contamination":{"+1":{"value":0.00000000,"err":0.00000000},"-1":{"value":0.00000000,"err":0.00000000}},"N_corrected":{"helicity":{"+1":{"value":57.00000000,"err":7.54983444},"-1":{"value":46.00000000,"err":6.78232998}},"total":{"value":103.00000000,"err":10.14889157}}},
    "(0,0,2,3)": {"N_data":{"helicity":{"+1":52,"-1":41},"total":93},"contamination":{"+1":{"value":0.03557312,"err":0.02325396},"-1":{"value":0.10358056,"err":0.05626935}},"N_corrected":{"helicity":{"+1":{"value":50.15019776,"err":7.05892185},"-1":{"value":36.75319704,"err":6.18617243}},"total":{"value":86.90339480,"err":9.38600592}}},
    "(0,0,2,4)": {"N_data":{"helicity":{"+1":142,"-1":123},"total":265},"contamination":{"+1":{"value":0.10493524,"err":0.02861631},"-1":{"value":0.08749254,"err":0.02631165}},"N_corrected":{"helicity":{"+1":{"value":127.09919592,"err":11.41377123},"-1":{"value":112.23841758,"err":10.62507620}},"total":{"value":239.33761350,"err":15.59379422}}},
    "(0,0,2,5)": {"N_data":{"helicity":{"+1":151,"-1":121},"total":272},"contamination":{"+1":{"value":0.07190083,"err":0.02060698},"-1":{"value":0.05454545,"err":0.01739696}},"N_corrected":{"helicity":{"+1":{"value":140.14297467,"err":11.82154681},"-1":{"value":114.40000055,"err":10.61089824}},"total":{"value":254.54297522,"err":15.88521735}}},
    "(0,0,2,6)": {"N_data":{"helicity":{"+1":94,"-1":94},"total":188},"contamination":{"+1":{"value":0.02865516,"err":0.01268671},"-1":{"value":0.06277743,"err":0.02230288}},"N_corrected":{"helicity":{"+1":{"value":91.30641496,"err":9.49274419},"-1":{"value":88.09892158,"err":9.32542155}},"total":{"value":179.40533654,"err":13.30697860}}},
    "(0,0,2,7)": {"N_data":{"helicity":{"+1":32,"-1":48},"total":80},"contamination":{"+1":{"value":0.09259259,"err":0.05024327},"-1":{"value":0.07834758,"err":0.03872684}},"N_corrected":{"helicity":{"+1":{"value":29.03703712,"err":5.37897705},"-1":{"value":44.23931616,"err":6.65046905}},"total":{"value":73.27635328,"err":8.55348658}}},
    "(0,0,2,8)": {"N_data":{"helicity":{"+1":24,"-1":24},"total":48},"contamination":{"+1":{"value":0.10869565,"err":0.08013150},"-1":{"value":0.11363636,"err":0.08392602}},"N_corrected":{"helicity":{"+1":{"value":21.39130440,"err":4.77123586},"-1":{"value":21.27272736,"err":4.78669743}},"total":{"value":42.66403176,"err":6.75848829}}},
    "(0,0,2,9)": {"N_data":{"helicity":{"+1":49,"-1":74},"total":123},"contamination":{"+1":{"value":0.00000000,"err":0.00000000},"-1":{"value":0.00000000,"err":0.00000000}},"N_corrected":{"helicity":{"+1":{"value":49.00000000,"err":7.00000000},"-1":{"value":74.00000000,"err":8.60232527}},"total":{"value":123.00000000,"err":11.09053651}}},
    "(0,0,2,10)": {"N_data":{"helicity":{"+1":75,"-1":95},"total":170},"contamination":{"+1":{"value":0.03825137,"err":0.04673221},"-1":{"value":0.00803213,"err":0.01091123}},"N_corrected":{"helicity":{"+1":{"value":72.13114725,"err":9.03639676},"-1":{"value":94.23694765,"err":9.72391357}},"total":{"value":166.36809490,"err":13.27444769}}},
    "(0,0,2,11)": {"N_data":{"helicity":{"+1":28,"-1":43},"total":71},"contamination":{"+1":{"value":0.00000000,"err":0.00000000},"-1":{"value":0.00000000,"err":0.00000000}},"N_corrected":{"helicity":{"+1":{"value":28.00000000,"err":5.29150262},"-1":{"value":43.00000000,"err":6.55743852}},"total":{"value":71.00000000,"err":8.42614977}}},
    "(0,0,3,0)": {"N_data":{"helicity":{"+1":17,"-1":14},"total":31},"contamination":{"+1":{"value":0.00000000,"err":0.00000000},"-1":{"value":0.00000000,"err":0.00000000}},"N_corrected":{"helicity":{"+1":{"value":17.00000000,"err":4.12310563},"-1":{"value":14.00000000,"err":3.74165739}},"total":{"value":31.00000000,"err":5.56776436}}},
    "(0,0,3,1)": {"N_data":{"helicity":{"+1":65,"-1":40},"total":105},"contamination":{"+1":{"value":0.02515723,"err":0.02145311},"-1":{"value":0.06410256,"err":0.05446058}},"N_corrected":{"helicity":{"+1":{"value":63.36478005,"err":7.98217981},"-1":{"value":37.43589760,"err":6.30727266}},"total":{"value":100.80067765,"err":10.17334177}}},
    "(0,0,3,2)": {"N_data":{"helicity":{"+1":41,"-1":33},"total":74},"contamination":{"+1":{"value":0.00000000,"err":0.00000000},"-1":{"value":0.00000000,"err":0.00000000}},"N_corrected":{"helicity":{"+1":{"value":41.00000000,"err":6.40312424},"-1":{"value":33.00000000,"err":5.74456265}},"total":{"value":74.00000000,"err":8.60232527}}},
    "(0,0,3,3)": {"N_data":{"helicity":{"+1":70,"-1":47},"total":117},"contamination":{"+1":{"value":0.07736944,"err":0.03498420},"-1":{"value":0.10101010,"err":0.04639428}},"N_corrected":{"helicity":{"+1":{"value":64.58413920,"err":8.09841851},"-1":{"value":42.25252530,"err":6.53753085}},"total":{"value":106.83666450,"err":10.40786683}}},
    "(0,0,3,4)": {"N_data":{"helicity":{"+1":224,"-1":155},"total":379},"contamination":{"+1":{"value":0.08658920,"err":0.01820458},"-1":{"value":0.12013305,"err":0.02627191}},"N_corrected":{"helicity":{"+1":{"value":204.60401920,"err":14.26590989},"-1":{"value":136.37937725,"err":11.68666249}},"total":{"value":340.98339645,"err":18.44164486}}},
    "(0,0,3,5)": {"N_data":{"helicity":{"+1":209,"-1":214},"total":423},"contamination":{"+1":{"value":0.05448029,"err":0.01255860},"-1":{"value":0.05555556,"err":0.01262229}},"N_corrected":{"helicity":{"+1":{"value":197.61361939,"err":13.91893929},"-1":{"value":202.11111016,"err":14.07760754}},"total":{"value":399.72472955,"err":19.79686604}}},
    "(0,0,3,6)": {"N_data":{"helicity":{"+1":104,"-1":92},"total":196},"contamination":{"+1":{"value":0.11263736,"err":0.03229066},"-1":{"value":0.14834267,"err":0.04177123}},"N_corrected":{"helicity":{"+1":{"value":92.28571456,"err":9.65238804},"-1":{"value":78.35247436,"err":9.02761032}},"total":{"value":170.63818892,"err":13.21613949}}},
    "(0,0,3,7)": {"N_data":{"helicity":{"+1":35,"-1":63},"total":98},"contamination":{"+1":{"value":0.18927126,"err":0.07848262},"-1":{"value":0.11670481,"err":0.04484292}},"N_corrected":{"helicity":{"+1":{"value":28.37550590,"err":5.52722825},"-1":{"value":55.64759697,"err":7.55873449}},"total":{"value":84.02310287,"err":9.36401192}}},
    "(0,0,3,8)": {"N_data":{"helicity":{"+1":57,"-1":49},"total":106},"contamination":{"+1":{"value":0.07173913,"err":0.03202772},"-1":{"value":0.11956522,"err":0.05233565}},"N_corrected":{"helicity":{"+1":{"value":52.91086959,"err":7.24208757},"-1":{"value":43.14130422,"err":6.67528968}},"total":{"value":96.05217381,"err":9.84922965}}},
    "(0,0,3,9)": {"N_data":{"helicity":{"+1":50,"-1":70},"total":120},"contamination":{"+1":{"value":0.04395604,"err":0.03161779},"-1":{"value":0.05023548,"err":0.03452502}},"N_corrected":{"helicity":{"+1":{"value":47.80219800,"err":6.94263742},"-1":{"value":66.48351640,"err":8.30568314}},"total":{"value":114.28571440,"err":10.82518298}}},
    "(0,0,3,10)": {"N_data":{"helicity":{"+1":53,"-1":63},"total":116},"contamination":{"+1":{"value":0.03846154,"err":0.03477118},"-1":{"value":0.02717391,"err":0.02500501}},"N_corrected":{"helicity":{"+1":{"value":50.96153838,"err":7.23862268},"-1":{"value":61.28804367,"err":7.88062353}},"total":{"value":112.24958205,"err":10.70055539}}},
    "(0,0,3,11)": {"N_data":{"helicity":{"+1":18,"-1":14},"total":32},"contamination":{"+1":{"value":0.01333333,"err":0.02007394},"-1":{"value":0.00000000,"err":0.00000000}},"N_corrected":{"helicity":{"+1":{"value":17.76000006,"err":4.20163779},"-1":{"value":14.00000000,"err":3.74165739}},"total":{"value":31.76000006,"err":5.62616745}}},
    "(0,0,4,1)": {"N_data":{"helicity":{"+1":24,"-1":25},"total":49},"contamination":{"+1":{"value":0.02222222,"err":0.02925692},"-1":{"value":0.03968254,"err":0.04719196}},"N_corrected":{"helicity":{"+1":{"value":23.46666672,"err":4.84130380},"-1":{"value":24.00793650,"err":4.94440758}},"total":{"value":47.47460322,"err":6.91992693}}},
    "(0,0,4,2)": {"N_data":{"helicity":{"+1":14,"-1":4},"total":18},"contamination":{"+1":{"value":0.00000000,"err":0.00000000},"-1":{"value":0.00000000,"err":0.00000000}},"N_corrected":{"helicity":{"+1":{"value":14.00000000,"err":3.74165739},"-1":{"value":4.00000000,"err":2.00000000}},"total":{"value":18.00000000,"err":4.24264069}}},
    "(0,0,4,3)": {"N_data":{"helicity":{"+1":55,"-1":26},"total":81},"contamination":{"+1":{"value":0.08750000,"err":0.03591258},"-1":{"value":0.16666667,"err":0.07672186}},"N_corrected":{"helicity":{"+1":{"value":50.18750000,"err":7.04964374},"-1":{"value":21.66666658,"err":4.69410867}},"total":{"value":71.85416658,"err":8.46948245}}},
    "(0,0,4,4)": {"N_data":{"helicity":{"+1":165,"-1":149},"total":314},"contamination":{"+1":{"value":0.13071437,"err":0.02704635},"-1":{"value":0.11204089,"err":0.02408589}},"N_corrected":{"helicity":{"+1":{"value":143.43212895,"err":12.02492052},"-1":{"value":132.30590739,"err":11.41760501}},"total":{"value":275.73803634,"err":16.58193045}}},
    "(0,0,4,5)": {"N_data":{"helicity":{"+1":179,"-1":165},"total":344},"contamination":{"+1":{"value":0.09501425,"err":0.01989543},"-1":{"value":0.09130361,"err":0.02053336}},"N_corrected":{"helicity":{"+1":{"value":161.99244925,"err":12.62076018},"-1":{"value":149.93490435,"err":12.15417131}},"total":{"value":311.92735360,"err":17.52162857}}},
    "(0,0,4,6)": {"N_data":{"helicity":{"+1":61,"-1":75},"total":136},"contamination":{"+1":{"value":0.12151865,"err":0.03755981},"-1":{"value":0.08806528,"err":0.02785346}},"N_corrected":{"helicity":{"+1":{"value":53.58736235,"err":7.23359243},"-1":{"value":68.39510400,"err":8.16920013}},"total":{"value":121.98246635,"err":10.91149349}}},
    "(0,0,4,7)": {"N_data":{"helicity":{"+1":33,"-1":46},"total":79},"contamination":{"+1":{"value":0.34567901,"err":0.11566928},"-1":{"value":0.17142857,"err":0.06020601}},"N_corrected":{"helicity":{"+1":{"value":21.59259267,"err":5.35711060},"-1":{"value":38.11428578,"err":6.26501462}},"total":{"value":59.70687845,"err":8.24312090}}},
    "(0,0,4,8)": {"N_data":{"helicity":{"+1":63,"-1":109},"total":172},"contamination":{"+1":{"value":0.20158103,"err":0.05826869},"-1":{"value":0.13586414,"err":0.03599037}},"N_corrected":{"helicity":{"+1":{"value":50.30039511,"err":7.32369430},"-1":{"value":94.19080874,"err":9.83784503}},"total":{"value":144.49120385,"err":12.26457064}}},
    "(0,0,4,9)": {"N_data":{"helicity":{"+1":24,"-1":46},"total":70},"contamination":{"+1":{"value":0.14210526,"err":0.08764674},"-1":{"value":0.02432432,"err":0.01847732}},"N_corrected":{"helicity":{"+1":{"value":20.58947376,"err":4.69983032},"-1":{"value":44.88108128,"err":6.67171686}},"total":{"value":65.47055504,"err":8.16089523}}},
    "(0,0,4,10)": {"N_data":{"helicity":{"+1":28,"-1":30},"total":58},"contamination":{"+1":{"value":0.00000000,"err":0.00000000},"-1":{"value":0.00000000,"err":0.00000000}},"N_corrected":{"helicity":{"+1":{"value":28.00000000,"err":5.29150262},"-1":{"value":30.00000000,"err":5.47722558}},"total":{"value":58.00000000,"err":7.61577311}}},
    "(0,0,4,11)": {"N_data":{"helicity":{"+1":0,"-1":2},"total":2},"contamination":{"+1":{"value":0.00000000,"err":0.00000000},"-1":{"value":0.00000000,"err":0.00000000}},"N_corrected":{"helicity":{"+1":{"value":0.00000000,"err":0.00000000},"-1":{"value":2.00000000,"err":1.41421356}},"total":{"value":2.00000000,"err":1.41421356}}},
    "(0,0,5,0)": {"N_data":{"helicity":{"+1":2,"-1":1},"total":3},"contamination":{"+1":{"value":0.00000000,"err":0.00000000},"-1":{"value":0.00000000,"err":0.00000000}},"N_corrected":{"helicity":{"+1":{"value":2.00000000,"err":1.41421356},"-1":{"value":1.00000000,"err":1.00000000}},"total":{"value":3.00000000,"err":1.73205081}}},
    "(0,0,5,1)": {"N_data":{"helicity":{"+1":27,"-1":13},"total":40},"contamination":{"+1":{"value":0.00000000,"err":0.00000000},"-1":{"value":0.00000000,"err":0.00000000}},"N_corrected":{"helicity":{"+1":{"value":27.00000000,"err":5.19615242},"-1":{"value":13.00000000,"err":3.60555128}},"total":{"value":40.00000000,"err":6.32455532}}},
    "(0,0,5,2)": {"N_data":{"helicity":{"+1":15,"-1":18},"total":33},"contamination":{"+1":{"value":0.08391608,"err":0.07665810},"-1":{"value":0.09090909,"err":0.08336463}},"N_corrected":{"helicity":{"+1":{"value":13.74125880,"err":3.72965825},"-1":{"value":16.36363638,"err":4.13856538}},"total":{"value":30.10489518,"err":5.57118247}}},
    "(0,0,5,3)": {"N_data":{"helicity":{"+1":25,"-1":21},"total":46},"contamination":{"+1":{"value":0.17647059,"err":0.08233360},"-1":{"value":0.07958478,"err":0.05369720}},"N_corrected":{"helicity":{"+1":{"value":20.58823525,"err":4.60345314},"-1":{"value":19.32871962,"err":4.36600760}},"total":{"value":39.91695487,"err":6.34458850}}},
    "(0,0,5,4)": {"N_data":{"helicity":{"+1":92,"-1":106},"total":198},"contamination":{"+1":{"value":0.16299796,"err":0.03836971},"-1":{"value":0.10518168,"err":0.02575948}},"N_corrected":{"helicity":{"+1":{"value":77.00418768,"err":8.77004311},"-1":{"value":94.85074192,"err":9.60884168}},"total":{"value":171.85492960,"err":13.00936181}}},
    "(0,0,5,5)": {"N_data":{"helicity":{"+1":110,"-1":103},"total":213},"contamination":{"+1":{"value":0.14407895,"err":0.03202494},"-1":{"value":0.11300310,"err":0.02645279}},"N_corrected":{"helicity":{"+1":{"value":94.15131550,"err":9.64343372},"-1":{"value":91.36068070,"err":9.40533303}},"total":{"value":185.51199620,"err":13.47056433}}},
    "(0,0,5,6)": {"N_data":{"helicity":{"+1":43,"-1":39},"total":82},"contamination":{"+1":{"value":0.12141816,"err":0.04403397},"-1":{"value":0.13497653,"err":0.04750212}},"N_corrected":{"helicity":{"+1":{"value":37.77901912,"err":6.06441700},"-1":{"value":33.73591533,"err":5.71090371}},"total":{"value":71.51493445,"err":8.33016055}}},
    "(0,0,5,7)": {"N_data":{"helicity":{"+1":29,"-1":29},"total":58},"contamination":{"+1":{"value":0.19115044,"err":0.07551574},"-1":{"value":0.19308126,"err":0.07376886}},"N_corrected":{"helicity":{"+1":{"value":23.45663724,"err":4.87532563},"-1":{"value":23.40064346,"err":4.84345013}},"total":{"value":46.85728070,"err":6.87224921}}},
    "(0,0,5,8)": {"N_data":{"helicity":{"+1":80,"-1":136},"total":216},"contamination":{"+1":{"value":0.25918704,"err":0.06138748},"-1":{"value":0.09332261,"err":0.02340471}},"N_corrected":{"helicity":{"+1":{"value":59.26503680,"err":8.24755798},"-1":{"value":123.30812504,"err":11.04230212}},"total":{"value":182.57316184,"err":13.78240359}}},
    "(0,0,5,9)": {"N_data":{"helicity":{"+1":46,"-1":67},"total":113},"contamination":{"+1":{"value":0.31231231,"err":0.12326107},"-1":{"value":0.17209046,"err":0.07079385}},"N_corrected":{"helicity":{"+1":{"value":31.63363374,"err":7.34187068},"-1":{"value":55.46993918,"err":8.27175458}},"total":{"value":87.10357292,"err":11.06006279}}},
    "(0,0,5,10)": {"N_data":{"helicity":{"+1":24,"-1":36},"total":60},"contamination":{"+1":{"value":0.00000000,"err":0.00000000},"-1":{"value":0.00000000,"err":0.00000000}},"N_corrected":{"helicity":{"+1":{"value":24.00000000,"err":4.89897949},"-1":{"value":36.00000000,"err":6.00000000}},"total":{"value":60.00000000,"err":7.74596669}}},
    "(0,0,5,11)": {"N_data":{"helicity":{"+1":1,"-1":2},"total":3},"contamination":{"+1":{"value":0.00000000,"err":0.00000000},"-1":{"value":0.00000000,"err":0.00000000}},"N_corrected":{"helicity":{"+1":{"value":1.00000000,"err":1.00000000},"-1":{"value":2.00000000,"err":1.41421356}},"total":{"value":3.00000000,"err":1.73205081}}},
    "(0,1,0,0)": {"N_data":{"helicity":{"+1":0,"-1":1},"total":1},"contamination":{"+1":{"value":0.00000000,"err":0.00000000},"-1":{"value":0.00000000,"err":0.00000000}},"N_corrected":{"helicity":{"+1":{"value":0.00000000,"err":0.00000000},"-1":{"value":1.00000000,"err":1.00000000}},"total":{"value":1.00000000,"err":1.00000000}}},
    "(0,1,0,1)": {"N_data":{"helicity":{"+1":2,"-1":1},"total":3},"contamination":{"+1":{"value":0.00000000,"err":0.00000000},"-1":{"value":0.00000000,"err":0.00000000}},"N_corrected":{"helicity":{"+1":{"value":2.00000000,"err":1.41421356},"-1":{"value":1.00000000,"err":1.00000000}},"total":{"value":3.00000000,"err":1.73205081}}},
    "(0,1,0,2)": {"N_data":{"helicity":{"+1":2,"-1":0},"total":2},"contamination":{"+1":{"value":0.00000000,"err":0.00000000},"-1":{"value":0.00000000,"err":0.00000000}},"N_corrected":{"helicity":{"+1":{"value":2.00000000,"err":1.41421356},"-1":{"value":0.00000000,"err":0.00000000}},"total":{"value":2.00000000,"err":1.41421356}}},
    "(0,1,0,4)": {"N_data":{"helicity":{"+1":3,"-1":5},"total":8},"contamination":{"+1":{"value":0.00000000,"err":0.00000000},"-1":{"value":0.00000000,"err":0.00000000}},"N_corrected":{"helicity":{"+1":{"value":3.00000000,"err":1.73205081},"-1":{"value":5.00000000,"err":2.23606798}},"total":{"value":8.00000000,"err":2.82842712}}},
    "(0,1,0,5)": {"N_data":{"helicity":{"+1":1,"-1":1},"total":2},"contamination":{"+1":{"value":0.00000000,"err":0.00000000},"-1":{"value":0.00000000,"err":0.00000000}},"N_corrected":{"helicity":{"+1":{"value":1.00000000,"err":1.00000000},"-1":{"value":1.00000000,"err":1.00000000}},"total":{"value":2.00000000,"err":1.41421356}}},
    "(0,1,0,6)": {"N_data":{"helicity":{"+1":5,"-1":0},"total":5},"contamination":{"+1":{"value":0.00000000,"err":0.00000000},"-1":{"value":0.00000000,"err":0.00000000}},"N_corrected":{"helicity":{"+1":{"value":5.00000000,"err":2.23606798},"-1":{"value":0.00000000,"err":0.00000000}},"total":{"value":5.00000000,"err":2.23606798}}},
    "(0,1,0,7)": {"N_data":{"helicity":{"+1":0,"-1":2},"total":2},"contamination":{"+1":{"value":0.00000000,"err":0.00000000},"-1":{"value":0.00000000,"err":0.00000000}},"N_corrected":{"helicity":{"+1":{"value":0.00000000,"err":0.00000000},"-1":{"value":2.00000000,"err":1.41421356}},"total":{"value":2.00000000,"err":1.41421356}}},
    "(0,1,0,8)": {"N_data":{"helicity":{"+1":0,"-1":1},"total":1},"contamination":{"+1":{"value":0.00000000,"err":0.00000000},"-1":{"value":0.00000000,"err":0.00000000}},"N_corrected":{"helicity":{"+1":{"value":0.00000000,"err":0.00000000},"-1":{"value":1.00000000,"err":1.00000000}},"total":{"value":1.00000000,"err":1.00000000}}},
    "(0,1,0,9)": {"N_data":{"helicity":{"+1":2,"-1":4},"total":6},"contamination":{"+1":{"value":0.00000000,"err":0.00000000},"-1":{"value":0.00000000,"err":0.00000000}},"N_corrected":{"helicity":{"+1":{"value":2.00000000,"err":1.41421356},"-1":{"value":4.00000000,"err":2.00000000}},"total":{"value":6.00000000,"err":2.44948974}}},
    "(0,1,0,10)": {"N_data":{"helicity":{"+1":1,"-1":2},"total":3},"contamination":{"+1":{"value":0.00000000,"err":0.00000000},"-1":{"value":0.00000000,"err":0.00000000}},"N_corrected":{"helicity":{"+1":{"value":1.00000000,"err":1.00000000},"-1":{"value":2.00000000,"err":1.41421356}},"total":{"value":3.00000000,"err":1.73205081}}},
    "(0,1,1,0)": {"N_data":{"helicity":{"+1":4,"-1":1},"total":5},"contamination":{"+1":{"value":0.00000000,"err":0.00000000},"-1":{"value":0.00000000,"err":0.00000000}},"N_corrected":{"helicity":{"+1":{"value":4.00000000,"err":2.00000000},"-1":{"value":1.00000000,"err":1.00000000}},"total":{"value":5.00000000,"err":2.23606798}}},
    "(0,1,1,1)": {"N_data":{"helicity":{"+1":13,"-1":6},"total":19},"contamination":{"+1":{"value":0.00000000,"err":0.00000000},"-1":{"value":0.00000000,"err":0.00000000}},"N_corrected":{"helicity":{"+1":{"value":13.00000000,"err":3.60555128},"-1":{"value":6.00000000,"err":2.44948974}},"total":{"value":19.00000000,"err":4.35889894}}},
    "(0,1,1,2)": {"N_data":{"helicity":{"+1":20,"-1":8},"total":28},"contamination":{"+1":{"value":0.00000000,"err":0.00000000},"-1":{"value":0.00000000,"err":0.00000000}},"N_corrected":{"helicity":{"+1":{"value":20.00000000,"err":4.47213595},"-1":{"value":8.00000000,"err":2.82842712}},"total":{"value":28.00000000,"err":5.29150262}}},
    "(0,1,1,3)": {"N_data":{"helicity":{"+1":4,"-1":5},"total":9},"contamination":{"+1":{"value":0.00000000,"err":0.00000000},"-1":{"value":0.00000000,"err":0.00000000}},"N_corrected":{"helicity":{"+1":{"value":4.00000000,"err":2.00000000},"-1":{"value":5.00000000,"err":2.23606798}},"total":{"value":9.00000000,"err":3.00000000}}},
    "(0,1,1,4)": {"N_data":{"helicity":{"+1":21,"-1":14},"total":35},"contamination":{"+1":{"value":0.06250000,"err":0.08975879},"-1":{"value":0.08333333,"err":0.12028131}},"N_corrected":{"helicity":{"+1":{"value":19.68750000,"err":4.69148267},"-1":{"value":12.83333338,"err":3.82093409}},"total":{"value":32.52083338,"err":6.05058237}}},
    "(0,1,1,5)": {"N_data":{"helicity":{"+1":18,"-1":14},"total":32},"contamination":{"+1":{"value":0.00000000,"err":0.00000000},"-1":{"value":0.00000000,"err":0.00000000}},"N_corrected":{"helicity":{"+1":{"value":18.00000000,"err":4.24264069},"-1":{"value":14.00000000,"err":3.74165739}},"total":{"value":32.00000000,"err":5.65685425}}},
    "(0,1,1,6)": {"N_data":{"helicity":{"+1":17,"-1":7},"total":24},"contamination":{"+1":{"value":0.02380952,"err":0.03692213},"-1":{"value":0.04761905,"err":0.07493293}},"N_corrected":{"helicity":{"+1":{"value":16.59523816,"err":4.07358453},"-1":{"value":6.66666665,"err":2.57377905}},"total":{"value":23.26190481,"err":4.81855056}}},
    "(0,1,1,7)": {"N_data":{"helicity":{"+1":4,"-1":3},"total":7},"contamination":{"+1":{"value":0.00000000,"err":0.00000000},"-1":{"value":0.00000000,"err":0.00000000}},"N_corrected":{"helicity":{"+1":{"value":4.00000000,"err":2.00000000},"-1":{"value":3.00000000,"err":1.73205081}},"total":{"value":7.00000000,"err":2.64575131}}},
    "(0,1,1,8)": {"N_data":{"helicity":{"+1":5,"-1":4},"total":9},"contamination":{"+1":{"value":0.00000000,"err":0.00000000},"-1":{"value":0.00000000,"err":0.00000000}},"N_corrected":{"helicity":{"+1":{"value":5.00000000,"err":2.23606798},"-1":{"value":4.00000000,"err":2.00000000}},"total":{"value":9.00000000,"err":3.00000000}}},
    "(0,1,1,9)": {"N_data":{"helicity":{"+1":13,"-1":22},"total":35},"contamination":{"+1":{"value":0.00000000,"err":0.00000000},"-1":{"value":0.00000000,"err":0.00000000}},"N_corrected":{"helicity":{"+1":{"value":13.00000000,"err":3.60555128},"-1":{"value":22.00000000,"err":4.69041576}},"total":{"value":35.00000000,"err":5.91607978}}},

etc etc.

yes I'd like you to try to write a fully fixed set of files everything we've discussed. Note that I actually just noticed I also get this (when I proceed to the bsa calculation phase which I hadn't looked at before):

Starting DVCS analysis...
Output directories ready.
Loaded binning scheme: 147 bins
[Loaded] sp18_inb, sp18_out, fa18_inb_supp, fa18_inb, fa18_out, sp19_inb
[Loaded] sp18_inb_gen, sp18_out_gen, fa18_inb_gen, fa18_out_gen, sp19_inb_gen
[Loaded] sp18_inb_rec, sp18_out_rec, fa18_inb_rec, fa18_out_rec, sp19_inb_rec
[Loaded] sp18_inb_gen_rad, sp18_out_gen_rad, fa18_inb_gen_rad, fa18_out_gen_rad, sp19_inb_gen_rad
[Loaded] sp18_inb_rec_rad, sp18_out_rec_rad, fa18_inb_rec_rad, fa18_out_rec_rad, sp19_inb_rec_rad
[Loaded] sp18_inb_eppi0, sp18_out_eppi0, fa18_inb_supp_eppi0, fa18_inb_eppi0, fa18_out_eppi0, sp19_inb_eppi0
[Loaded] sp18_inb_gen_mc, sp18_out_gen_mc, fa18_inb_gen_mc, fa18_out_gen_mc, sp19_inb_gen_mc
[Loaded] sp18_inb_rec_mc, sp18_out_rec_mc, fa18_inb_rec_mc, fa18_out_rec_mc, sp19_inb_rec_mc
[Loaded] sp18_inb_bkg, sp18_out_bkg, fa18_inb_bkg, fa18_out_bkg, sp19_inb_bkg
All trees loaded successfully.
[pi0_corr] Wrote output/jsons/pi0_corrected_counts_DVCS_Fa18_inb.json
[pi0_corr] Wrote output/jsons/pi0_corrected_counts_DVCS_Fa18_out.json
[pi0_corr] Wrote output/jsons/pi0_corrected_counts_DVCS_Sp19_inb.json
[pi0_corr] Wrote output/jsons/pi0_corrected_counts_DVCS_Sp18_out.json
[pi0_corr] Wrote output/jsons/pi0_corrected_counts_DVCS_Sp18_inb.json
[pi0_corr] Wrote output/jsons/pi0_corrected_counts_DVCS_Fa18_inb_supp.json
[pi0_corr] Wrote combined output/jsons/pi0_corrected_counts_all_periods.json
π0-corrected counts stage finished.
[bsa][WARN] total_counts has no group 'fa18_inb'. Skipping DVCS_Fa18_inb
[bsa][WARN] total_counts has no group 'fa18_out'. Skipping DVCS_Fa18_out
[bsa][WARN] total_counts has no group 'sp19_inb'. Skipping DVCS_Sp19_inb
[bsa][WARN] total_counts has no group 'sp18_out'. Skipping DVCS_Sp18_out
[bsa][WARN] total_counts has no group 'sp18_inb'. Skipping DVCS_Sp18_inb
[bsa][WARN] total_counts has no group 'fa18_inb_supp'. Skipping DVCS_Fa18_inb_supp
Info in <TCanvas::Print>: png file output/bsa_plots/10.6_combined/plot_bsa_RGA_10.6_combined_xB_0.png has been created
Info in <TCanvas::Print>: png file output/bsa_plots/10.6_combined/plot_bsa_RGA_10.6_combined_xB_1.png has been created
Info in <TCanvas::Print>: png file output/bsa_plots/10.6_combined/plot_bsa_RGA_10.6_combined_xB_2.png has been created
Info in <TCanvas::Print>: png file output/bsa_plots/10.6_combined/plot_bsa_RGA_10.6_combined_xB_3.png has been created
Info in <TCanvas::Print>: png file output/bsa_plots/10.6_combined/plot_bsa_RGA_10.6_combined_xB_4.png has been created
Info in <TCanvas::Print>: png file output/bsa_plots/10.6_combined/plot_bsa_RGA_10.6_combined_xB_5.png has been created
Info in <TCanvas::Print>: png file output/bsa_plots/10.6_combined/plot_bsa_RGA_10.6_combined_xB_6.png has been created
Info in <TCanvas::Print>: png file output/bsa_plots/10.6_combined/plot_bsa_RGA_10.6_combined_xB_7.png has been created
All done.

i.e. the total counts might have the wrong names? anyway, here is my current main.cpp:

#include "make_dirs.h"
#include "load_trees.h"
#include "exclusivity_cuts.h"
#include "load_binning_scheme.h"
#include "bin_means.h"
#include "total_counts.h"
#include "pi0_contamination.h"
#include "pi0_corrected_counts.h"
#include "bsa.h"
#include "radiative_corrections.h"  
#include <filesystem>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include "acceptance.h"
#include "unfolding.h"
#include "bin_volume.h"
#include "uncorrected_cross_section.h"
#include "rad_corrected_cross_section.h"

int main(int argc, char* argv[]) {
    std::cout << "Starting DVCS analysis..." << std::endl;

    // Create necessary output directories
    makeOutputDirs();
    std::cout << "Output directories ready." << std::endl;

    // Root of output tree (used by several stages)
    const std::string output_root = "output";

    // --- Load binning scheme ---
    const std::string csv_file_path = "imports/integrated_bin_v2.csv";
    auto binning_scheme = load_binning_scheme(csv_file_path);
    std::cout << "Loaded binning scheme: " << binning_scheme.size() << " bins" << std::endl;

    // Containers for different tree categories
    std::map<std::string, TTree*> dataTrees;        // DVCS data
    std::map<std::string, TTree*> genMcTrees;       // DVCS generated MC (no-rad)
    std::map<std::string, TTree*> recMcTrees;       // DVCS reconstructed MC (no-rad)
    std::map<std::string, TTree*> eppi0DataTrees;   // eppi0 data
    std::map<std::string, TTree*> eppi0GenMcTrees;  // eppi0 generated MC
    std::map<std::string, TTree*> eppi0RecMcTrees;  // eppi0 reconstructed MC
    std::map<std::string, TTree*> radGenMcTrees;    // NEW: DVCS generated MC (radiative)
    std::map<std::string, TTree*> radRecMcTrees;    // NEW: DVCS reconstructed MC (radiative)

    // Load all trees from files
    loadTrees(dataTrees, genMcTrees, recMcTrees,
              eppi0DataTrees, eppi0GenMcTrees, eppi0RecMcTrees,
              radGenMcTrees, radRecMcTrees);
    std::cout << "All trees loaded successfully." << std::endl;

    // // Run exclusivity cut extraction (single-threaded for stability)
    // runAllExclusivityCuts(
    //     dataTrees, recMcTrees, eppi0DataTrees, eppi0RecMcTrees,
    //     "output/jsons", "output/exclusivity_plots", 1
    // );
    // std::cout << "Exclusivity-cut stage finished." << std::endl;

    // --------- Global bin-averaged kinematics ----------
    std::vector<std::string> dvcs_periods = {
        "DVCS_Fa18_inb", "DVCS_Fa18_out", "DVCS_Sp19_inb",
        "DVCS_Sp18_out", "DVCS_Sp18_inb", "DVCS_Fa18_inb_supp"
    };
    std::vector<std::string> topologies = {"(FD,FD)", "(CD,FD)", "(CD,FT)"};
    const std::string analysis_type = "dvcs";
    const std::string output_json_means = "output/jsons/bin_means_global.json";

    // calculate_bin_means(dvcs_periods, topologies, analysis_type, binning_scheme, output_json_means, 
    //     dataTrees);

    // --------- Total counts after exclusivity cuts (by helicity) ----------
    const std::string cuts_json_path   = "output/jsons/combined_cuts.json"; 
    // produced by exclusivity_cuts
    const std::string output_counts_js = "output/jsons/total_counts.json";

    // compute_total_counts(dvcs_periods, topologies, binning_scheme, dataTrees, cuts_json_path, 
    //     output_counts_js);

    // Helicity-resolved π0 contamination
    // NOTE: pass the OUTPUT ROOT ("output") so the implementation writes:
    //   - per-period JSONs to output/jsons/contamination/
    //   - combined JSON to output/jsons/
    //   - plots to output/contamination_plots/...
    // compute_pi0_contamination_helicity(
    //     dvcs_periods,
    //     topologies,
    //     binning_scheme,
    //     dataTrees,
    //     eppi0DataTrees,
    //     eppi0RecMcTrees,   // reco MC (keys "*_rec_mc")
    //     eppi0RecMcTrees,   // bkg MC  (keys "*_bkg")
    //     cuts_json_path,
    //     output_root        // <<< was "output/contamination"; must be the root "output"
    // );

    // --------- π0-corrected helicity counts (per φ) ----------
    const std::string total_counts_json = "output/jsons/total_counts.json";
    const std::string contamination_dir_counts = "output/jsons/contamination";
    compute_pi0_corrected_counts(
        dvcs_periods,
        binning_scheme,
        total_counts_json,        // from compute_total_counts()
        contamination_dir_counts, // from compute_pi0_contamination_helicity()
        output_root               // "output"
    );
    std::cout << "π0-corrected counts stage finished." << std::endl;

    // --------- Radiative corrections (generated Born/Rad, per bin & φ) ----------
    // Writes:
    //   - per-group JSONs: output/jsons/radiative_corrections_group_<energy>.json
    //   - per-period JSONs: output/jsons/radiative_corrections_<period>.json
    //   - all-groups file: output/jsons/radiative_corrections_all_groups.json
    //   - plots (ONLY per beam energy): output/radiative_correction_plots/{10.59,10.60,10.2}/...
    // compute_radiative_corrections(
    //     dvcs_periods,
    //     binning_scheme,
    //     genMcTrees,
    //     radGenMcTrees,
    //     output_root
    // );

    // Beam-Spin Asymmetry:
    // - reads total_counts.json and contamination JSONs
    // - writes per-period fits to output/jsons/BSA_fits/BSA_fits_<period>.json
    // - writes all-periods file to output/jsons/BSA_fits_all_periods.json
    // - writes 10.6 GeV combined to output/jsons/BSA_fits_combined_10p6.json
    // - plots to output/bsa_plots/<runTag>/...
    namespace fs = std::filesystem;
    const std::string contamination_dir_bsa =
        (fs::path(output_root) / "jsons" / "contamination").string();
    compute_and_plot_bsa_helicity(
        dvcs_periods,
        topologies,
        binning_scheme,
        dataTrees,            // DVCS trees (for beam_pol extraction)
        output_counts_js,     // total_counts.json path
        contamination_dir_bsa,// directory with contamination_<period>.json files
        output_root
    );

    //     // --------- Acceptance (reco MC with MC cuts / gen MC, per bin & φ) ----------
    // // Writes per-period JSONs: output/jsons/acceptance_<period>.json
    // // Plots to: output/acceptance/<runTag>/plot_acceptance_<period>_xB_<ix>.png
    // {
    //     std::vector<std::string> acc_periods = {
    //         "DVCS_Sp18_inb", "DVCS_Sp18_out",
    //         "DVCS_Fa18_inb", "DVCS_Fa18_out",
    //         "DVCS_Sp19_inb"
    //     }; // intentionally skipping DVCS_Fa18_inb_supp
    //     compute_and_plot_acceptance(
    //         acc_periods,
    //         topologies,
    //         binning_scheme,
    //         genMcTrees,
    //         recMcTrees,
    //         cuts_json_path,
    //         output_root
    //     );
    // }

    // // --------- Unfolding (counts / acceptance), helicity-resolved ----------
    // // Writes per-period JSONs: output/jsons/unfolded_<period>.json
    // // Plots to: output/unfolding/<runTag>/plot_unfolded_<period>_xB_<ix>.png
    // {
    //     std::vector<std::string> unf_periods = {
    //         "DVCS_Sp18_inb", "DVCS_Sp18_out",
    //         "DVCS_Fa18_inb", "DVCS_Fa18_out",
    //         "DVCS_Sp19_inb"
    //     }; // skip DVCS_Fa18_inb_supp on purpose
    //     const std::string total_counts_js = "output/jsons/total_counts.json";
    //     compute_and_plot_unfolding(
    //         unf_periods,
    //         binning_scheme,
    //         total_counts_js,
    //         output_root
    //     );
    // }

    // // --------- Bin Volume (generator-based φ coverage), per beam energy ----------
    // // Writes per-energy JSONs: output/jsons/bin_volume_<energy>.json
    // // Plots to: output/bin_volume/<energy>/plot_bin_volume_<energy>_xB_<ix>.png
    // compute_and_plot_bin_volume(
    //     binning_scheme,
    //     genMcTrees,
    //     output_root
    // );

    // compute_uncorrected_cross_sections(
    //     binning_scheme,
    //     "output/jsons",                    // bin volume JSON directory
    //     "output/jsons",               // unfolded counts per helicity
    //     "imports/integrated_luminosity",        // luminosity text files
    //     "output/uncorrected_cross_section"      // output dir
    // );

    // {
    //     // Multiply uncorrected dσ/dφ by Born/Rad per-φ correction
    //     const std::string unx_dir = "output/uncorrected_cross_section/jsons";
    //     const std::string rc_dir  = "output/jsons"; // radiative_corrections_group_<E>.json
    //     const std::string out_dir = "output/rad_corrected_cross_section";
    //     compute_rad_corrected_cross_sections(binning_scheme, unx_dir, rc_dir, out_dir);
    // }

    std::cout << "All done." << std::endl;
    return 0;
}

--

Here is pi0_contamination.cpp:

#include "pi0_contamination.h"

#include <TTree.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

bool COPY_CONTAM_TO_FA18_INB_SUPP = true;
bool ENABLE_PI0_CONTAMINATION_PLOTS = true;

// ---------------- helpers (shared style with total_counts.cpp) ----------------
static inline std::string toLower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c){ return std::tolower(c); });
    return s;
}
static std::string periodToRunTagKey(const std::string& period) {
    auto pos = period.find('_');
    if (pos == std::string::npos || pos + 1 >= period.size()) return toLower(period);
    return toLower(period.substr(pos + 1));
}
static inline std::string topoToKey(const std::string& topoStr) {
    if (topoStr == "(FD,FD)") return "FD_FD";
    if (topoStr == "(CD,FD)") return "CD_FD";
    if (topoStr == "(CD,FT)") return "CD_FT";
    return "FD_FD";
}
static inline bool passesTopology_simple(int detector1, int detector2, const std::string& topoStr) {
    if (topoStr == "(FD,FD)") return (detector1 == 1 && detector2 == 1);
    if (topoStr == "(CD,FD)") return (detector1 == 2 && detector2 == 1);
    if (topoStr == "(CD,FT)") return (detector1 == 2 && detector2 == 0);
    return false;
}
static inline bool applyKinematicCuts_simple(double t1, double open_angle_ep2, double pTmiss) {
    if (open_angle_ep2 <= 5.0) return false;
    if ((-t1) > 1.0)          return false;
    if (pTmiss > 0.20)        return false;
    return true;
}

// binning helpers
static constexpr int    N_PHI_BINS = 12;
static constexpr double TWO_PI     = 2.0 * M_PI;
static inline double wrapToTwoPi(double phi) {
    if (!std::isfinite(phi)) return 0.0;
    double w = std::fmod(phi, TWO_PI);
    if (w < 0.0) w += TWO_PI;
    if (w >= TWO_PI) w = std::nextafter(TWO_PI, 0.0);
    return w;
}
static inline int phiToBin(double phi) {
    double w = wrapToTwoPi(phi);
    double width = TWO_PI / static_cast<double>(N_PHI_BINS);
    int idx = static_cast<int>(std::floor(w / width));
    if (idx < 0) idx = 0;
    if (idx >= N_PHI_BINS) idx = N_PHI_BINS - 1;
    return idx;
}
static std::vector<std::pair<double,double>> uniqueRanges(const std::vector<Binning>& scheme, char which) {
    std::set<std::pair<double,double>> s;
    for (const auto& b : scheme) {
        if (which == 'x') s.emplace(b.xBmin, b.xBmax);
        else if (which == 'Q') s.emplace(b.Q2min, b.Q2max);
        else if (which == 't') s.emplace(b.tmin, b.tmax);
    }
    return std::vector<std::pair<double,double>>(s.begin(), s.end());
}
static int findBin(double v, const std::vector<std::pair<double,double>>& ranges) {
    for (int i = 0; i < static_cast<int>(ranges.size()); ++i)
        if (v >= ranges[i].first && v < ranges[i].second) return i;
    return -1;
}

// ---------------- exclusivity cuts loader (data-cuts only) ----------------
struct Stats { double mean = 0.0; double std = 0.0; };
using VarCutMap = std::map<std::string, Stats>;
using PeriodTopoCuts = std::map<std::string, VarCutMap>; // key: "DVCS_Fa18_inb_FD_FD"

static bool loadCombinedCuts(const std::string& path, PeriodTopoCuts& out) {
    std::ifstream ifs(path);
    if (!ifs) { std::cerr << "[pi0_contam][ERROR] Cannot open cuts JSON: " << path << std::endl; return false; }
    std::string s((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    ifs.close();

    size_t pos = 0;
    while (true) {
        size_t keyStart = s.find('"', pos);
        if (keyStart == std::string::npos) break;
        size_t keyEnd = s.find('"', keyStart + 1);
        if (keyEnd == std::string::npos) break;
        std::string key = s.substr(keyStart + 1, keyEnd - keyStart - 1);

        size_t dataPos = s.find("\"data\"", keyEnd);
        if (dataPos == std::string::npos) { pos = keyEnd + 1; continue; }
        size_t braceStart = s.find('{', dataPos);
        if (braceStart == std::string::npos) { pos = keyEnd + 1; continue; }

        int depth = 0; size_t i = braceStart;
        for (; i < s.size(); ++i) {
            if (s[i] == '{') depth++;
            else if (s[i] == '}') { depth--; if (depth == 0) { ++i; break; } }
        }
        if (depth != 0) { pos = keyEnd + 1; continue; }
        std::string dataObj = s.substr(braceStart, i - braceStart);

        VarCutMap cuts;
        size_t vpos = 0;
        while (true) {
            size_t vKeyS = dataObj.find('"', vpos);
            if (vKeyS == std::string::npos) break;
            size_t vKeyE = dataObj.find('"', vKeyS + 1);
            if (vKeyE == std::string::npos) break;
            std::string var = dataObj.substr(vKeyS + 1, vKeyE - vKeyS - 1);

            size_t meanPos = dataObj.find("\"mean\"", vKeyE);
            size_t stdPos  = dataObj.find("\"std\"",  vKeyE);
            if (meanPos == std::string::npos || stdPos == std::string::npos) { vpos = vKeyE + 1; continue; }

            auto readNum = [&](size_t from)->double {
                size_t colon = dataObj.find(':', from);
                if (colon == std::string::npos) return 0.0;
                size_t j = colon + 1;
                while (j < dataObj.size() && std::isspace(static_cast<unsigned char>(dataObj[j]))) ++j;
                size_t k = j;
                while (k < dataObj.size() && (std::isdigit(static_cast<unsigned char>(dataObj[k])) || dataObj[k]=='-' || dataObj[k]=='+' || dataObj[k]=='.' || dataObj[k]=='e' || dataObj[k]=='E')) ++k;
                try { return std::stod(dataObj.substr(j, k - j)); } catch (...) { return 0.0; }
            };

            double m = readNum(meanPos);
            double sd = readNum(stdPos);
            cuts[var] = Stats{m, sd};

            vpos = vKeyE + 1;
        }

        if (key.rfind("DVCS_", 0) == 0) out[key] = cuts; // keep only DVCS entries
        pos = keyEnd + 1;
    }
    return !out.empty();
}

static inline bool within3Sigma(double v, const Stats& s) {
    return (v >= s.mean - 3.0*s.std) && (v <= s.mean + 3.0*s.std);
}
static bool passes3SigmaCuts(const VarCutMap& cuts, const std::map<std::string,double>& values) {
    for (const auto& kv : cuts) {
        auto it = values.find(kv.first);
        if (it == values.end()) continue;
        if (!within3Sigma(it->second, kv.second)) return false;
    }
    return true;
}

// ---------------- branch binders ----------------
struct BranchBinderDVCS {
    int detector1=0, detector2=0; bool has_d1=false, has_d2=false;
    int helicity=0; bool has_helicity=false;
    double t1=0, open_angle_ep2=0, pTmiss=0; bool has_t1=false, has_oa=false, has_pT=false;
    double x=0, Q2=0, phi2=0, Delta_phi=0; bool has_x=false, has_Q2=false, has_phi2=false, has_Dp=false;
    // exclusivity vars
    double Emiss2=0, Mx2=0, Mx2_1=0, Mx2_2=0, theta_gamma_gamma=0, xF=0;
    bool has_Em=false, has_Mx2=false, has_Mx21=false, has_Mx22=false, has_tgg=false, has_xF=false;

    void bind(TTree* t) {
        if (!t) return;
        auto bindI=[&](const char*n,int*a,bool&f){ if(t->GetBranch(n)){ t->SetBranchAddress(n,a); f=true; } };
        auto bindD=[&](const char*n,double*a,bool&f){ if(t->GetBranch(n)){ t->SetBranchAddress(n,a); f=true; } };

        bindI("detector1",&detector1,has_d1);
        bindI("detector2",&detector2,has_d2);
        bindI("helicity",&helicity,has_helicity);

        bindD("t1",&t1,has_t1);
        bindD("open_angle_ep2",&open_angle_ep2,has_oa);
        bindD("pTmiss",&pTmiss,has_pT);

        bindD("x",&x,has_x);
        bindD("Q2",&Q2,has_Q2);
        bindD("phi2",&phi2,has_phi2);
        bindD("Delta_phi",&Delta_phi,has_Dp);

        bindD("Emiss2",&Emiss2,has_Em);
        bindD("Mx2",&Mx2,has_Mx2);
        bindD("Mx2_1",&Mx2_1,has_Mx21);
        bindD("Mx2_2",&Mx2_2,has_Mx22);
        bindD("theta_gamma_gamma",&theta_gamma_gamma,has_tgg);
        bindD("xF",&xF,has_xF);
    }
    bool readyCuts() const { return has_d1 && has_d2 && has_t1 && has_oa && has_pT && has_helicity; }
    bool readyBins() const { return has_x && has_Q2 && (has_phi2 || has_Dp); }
    double phi() const { return has_phi2 ? phi2 : (has_Dp ? Delta_phi : std::numeric_limits<double>::quiet_NaN()); }

    std::map<std::string,double> cutVals() const {
        std::map<std::string,double> m;
        if (has_Dp)  m["Delta_phi"]=Delta_phi;
        if (has_tgg) m["theta_gamma_gamma"]=theta_gamma_gamma;
        if (has_pT)  m["pTmiss"]=pTmiss;
        if (has_xF)  m["xF"]=xF;
        if (has_Em)  m["Emiss2"]=Emiss2;
        if (has_Mx2) m["Mx2"]=Mx2;
        if (has_Mx21) m["Mx2_1"]=Mx2_1;
        if (has_Mx22) m["Mx2_2"]=Mx2_2;
        return m;
    }
};
struct BranchBinderEPPI0Data { // has helicity
    int detector1=0, detector2=0; bool has_d1=false, has_d2=false;
    int helicity=0; bool has_helicity=false;
    double t1=0, open_angle_ep2=0, pTmiss=0; bool has_t1=false, has_oa=false, has_pT=false;
    double x=0, Q2=0, phi2=0, Delta_phi=0; bool has_x=false, has_Q2=false, has_phi2=false, has_Dp=false;
    // exclusivity vars use theta_pi0_pi0 for eppi0
    double Emiss2=0, Mx2=0, Mx2_1=0, Mx2_2=0, theta_pi0_pi0=0, xF=0;
    bool has_Em=false, has_Mx2=false, has_Mx21=false, has_Mx22=false, has_tpp=false, has_xF=false;

    void bind(TTree* t) {
        if (!t) return;
        auto bindI=[&](const char*n,int*a,bool&f){ if(t->GetBranch(n)){ t->SetBranchAddress(n,a); f=true; } };
        auto bindD=[&](const char*n,double*a,bool&f){ if(t->GetBranch(n)){ t->SetBranchAddress(n,a); f=true; } };

        bindI("detector1",&detector1,has_d1);
        bindI("detector2",&detector2,has_d2);
        bindI("helicity",&helicity,has_helicity);

        bindD("t1",&t1,has_t1);
        bindD("open_angle_ep2",&open_angle_ep2,has_oa);
        bindD("pTmiss",&pTmiss,has_pT);

        bindD("x",&x,has_x);
        bindD("Q2",&Q2,has_Q2);
        bindD("phi2",&phi2,has_phi2);
        bindD("Delta_phi",&Delta_phi,has_Dp);

        bindD("Emiss2",&Emiss2,has_Em);
        bindD("Mx2",&Mx2,has_Mx2);
        bindD("Mx2_1",&Mx2_1,has_Mx21);
        bindD("Mx2_2",&Mx2_2,has_Mx22);
        bindD("theta_pi0_pi0",&theta_pi0_pi0,has_tpp);
        bindD("xF",&xF,has_xF);
    }
    bool readyCuts() const { return has_d1 && has_d2 && has_t1 && has_oa && has_pT && has_helicity; }
    bool readyBins() const { return has_x && has_Q2 && (has_phi2 || has_Dp); }
    double phi() const { return has_phi2 ? phi2 : (has_Dp ? Delta_phi : std::numeric_limits<double>::quiet_NaN()); }

    std::map<std::string,double> cutVals() const {
        std::map<std::string,double> m;
        if (has_Dp)  m["Delta_phi"]=Delta_phi;
        if (has_tpp) m["theta_pi0_pi0"]=theta_pi0_pi0;
        if (has_pT)  m["pTmiss"]=pTmiss;
        if (has_xF)  m["xF"]=xF;
        if (has_Em)  m["Emiss2"]=Emiss2;
        if (has_Mx2) m["Mx2"]=Mx2;
        if (has_Mx21) m["Mx2_1"]=Mx2_1;
        if (has_Mx22) m["Mx2_2"]=Mx2_2;
        return m;
    }
};
struct BranchBinderEPPI0MC { // no helicity
    int detector1=0, detector2=0; bool has_d1=false, has_d2=false;
    double t1=0, open_angle_ep2=0, pTmiss=0; bool has_t1=false, has_oa=false, has_pT=false;
    double x=0, Q2=0, phi2=0, Delta_phi=0; bool has_x=false, has_Q2=false, has_phi2=false, has_Dp=false;
    // exclusivity
    double Emiss2=0, Mx2=0, Mx2_1=0, Mx2_2=0, theta_pi0_pi0=0, theta_gamma_gamma=0, xF=0;
    bool has_Em=false, has_Mx2=false, has_Mx21=false, has_Mx22=false, has_tpp=false, has_tgg=false, has_xF=false;

    void bind(TTree* t) {
        if (!t) return;
        auto bindI=[&](const char*n,int*a,bool&f){ if(t->GetBranch(n)){ t->SetBranchAddress(n,a); f=true; } };
        auto bindD=[&](const char*n,double*a,bool&f){ if(t->GetBranch(n)){ t->SetBranchAddress(n,a); f=true; } };

        bindI("detector1",&detector1,has_d1);
        bindI("detector2",&detector2,has_d2);

        bindD("t1",&t1,has_t1);
        bindD("open_angle_ep2",&open_angle_ep2,has_oa);
        bindD("pTmiss",&pTmiss,has_pT);

        bindD("x",&x,has_x);
        bindD("Q2",&Q2,has_Q2);
        bindD("phi2",&phi2,has_phi2);
        bindD("Delta_phi",&Delta_phi,has_Dp);

        bindD("Emiss2",&Emiss2,has_Em);
        bindD("Mx2",&Mx2,has_Mx2);
        bindD("Mx2_1",&Mx2_1,has_Mx21);
        bindD("Mx2_2",&Mx2_2,has_Mx22);
        bindD("theta_pi0_pi0",&theta_pi0_pi0,has_tpp);
        bindD("theta_gamma_gamma",&theta_gamma_gamma,has_tgg);
        bindD("xF",&xF,has_xF);
    }
    bool readyCuts() const { return has_d1 && has_d2 && has_t1 && has_oa && has_pT; }
    bool readyBins() const { return has_x && has_Q2 && (has_phi2 || has_Dp); }
    double phi() const { return has_phi2 ? phi2 : (has_Dp ? Delta_phi : std::numeric_limits<double>::quiet_NaN()); }

    std::map<std::string,double> cutValsForDVCS() const { // when mis-ID to DVCS hypothesis
        std::map<std::string,double> m;
        if (has_Dp)  m["Delta_phi"]=Delta_phi;
        if (has_tgg) m["theta_gamma_gamma"]=theta_gamma_gamma;
        if (has_pT)  m["pTmiss"]=pTmiss;
        if (has_xF)  m["xF"]=xF;
        if (has_Em)  m["Emiss2"]=Emiss2;
        if (has_Mx2) m["Mx2"]=Mx2;
        if (has_Mx21) m["Mx2_1"]=Mx2_1;
        if (has_Mx22) m["Mx2_2"]=Mx2_2;
        return m;
    }
    std::map<std::string,double> cutValsForEPPI0() const { // when genuine π0 selection
        std::map<std::string,double> m;
        if (has_Dp)  m["Delta_phi"]=Delta_phi;
        if (has_tpp) m["theta_pi0_pi0"]=theta_pi0_pi0;
        if (has_pT)  m["pTmiss"]=pTmiss;
        if (has_xF)  m["xF"]=xF;
        if (has_Em)  m["Emiss2"]=Emiss2;
        if (has_Mx2) m["Mx2"]=Mx2;
        if (has_Mx21) m["Mx2_1"]=Mx2_1;
        if (has_Mx22) m["Mx2_2"]=Mx2_2;
        return m;
    }
};

// ---------------- containers ----------------
using BinKey = std::tuple<int,int,int,int>; // (ix,iQ2,it,ip)
struct HelCounts { long long plus=0, minus=0; };
struct BinCounts {
    HelCounts N_data;       // DVCS data, by helicity
    HelCounts N_pi0_exp;    // eppi0 DATA, by helicity
    long long N_pi0_mc  = 0; // eppi0_bkg mis-ID MC (no helicity)
    long long N_pi0_reco= 0; // eppi0 reco MC (no helicity)
    // results (per helicity)
    double c_plus = 0.0, c_plus_err = 0.0;
    double c_minus= 0.0, c_minus_err= 0.0;
};

// ---------------- write JSON ----------------
static inline std::string keyStr(const BinKey& k) {
    int a,b,c,d; std::tie(a,b,c,d)=k;
    std::ostringstream os; os<<"("<<a<<","<<b<<","<<c<<","<<d<<")";
    return os.str();
}

static void writeContaminationJson(const std::string& path,
                                   const std::map<BinKey, BinCounts>& table,
                                   int nPhi,
                                   const std::vector<std::pair<double,double>>& xB_bins,
                                   const std::vector<std::pair<double,double>>& Q2_bins,
                                   const std::vector<std::pair<double,double>>& t_bins) {
    std::ofstream ofs(path);
    if (!ofs) { std::cerr << "[pi0_contam][ERROR] Cannot open " << path << std::endl; return; }
    ofs << std::fixed << std::setprecision(8);
    ofs << "{\n";
    ofs << "  \"binning_meta\": {\"phi_bins\": " << nPhi
        << ", \"xB_bins\": " << xB_bins.size()
        << ", \"Q2_bins\": " << Q2_bins.size()
        << ", \"t_bins\": "  << t_bins.size()  << "},\n";
    ofs << "  \"bins\": {\n";
    bool first=true;
    for (const auto& kv : table) {
        if (!first) ofs << ",\n";
        first=false;
        const auto& bc = kv.second;
        ofs << "    \"" << keyStr(kv.first) << "\": {"
            << "\"N_data\":{\"helicity\":{\"+1\":" << bc.N_data.plus
            << ",\"-1\":" << bc.N_data.minus << "},\"total\":" << (bc.N_data.plus+bc.N_data.minus) << "}"
            << ",\"N_pi0_exp\":{\"helicity\":{\"+1\":" << bc.N_pi0_exp.plus
            << ",\"-1\":" << bc.N_pi0_exp.minus << "},\"total\":" << (bc.N_pi0_exp.plus+bc.N_pi0_exp.minus) << "}"
            << ",\"N_pi0_mc\":"   << bc.N_pi0_mc
            << ",\"N_pi0_reco\":" << bc.N_pi0_reco
            << ",\"contamination\":{"
            << "\"+1\":{\"value\":" << bc.c_plus  << ",\"err\":" << bc.c_plus_err  << "},"
            << "\"-1\":{\"value\":" << bc.c_minus << ",\"err\":" << bc.c_minus_err << "}"
            << "}"
            << "}";
    }
    ofs << "\n  }\n}\n";
    ofs.close();
    std::cout << "[pi0_contam] Wrote " << path << std::endl;
}

// =====================================================
// ROOT plotting for helicity-resolved π0 contamination
// =====================================================
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TH1.h>
#include <TH1F.h>
#include <TLatex.h>
#include <TPad.h>

static constexpr int N_PHI_BINS_PLOT = 12;

// Phi centers in degrees (0..360), consistent with equal-width 12 bins
static std::vector<double> phiCentersDeg() {
    std::vector<double> v(N_PHI_BINS_PLOT);
    const double step = 360.0 / static_cast<double>(N_PHI_BINS_PLOT);
    for (int i = 0; i < N_PHI_BINS_PLOT; ++i) v[i] = (i + 0.5) * step;
    return v;
}

// Collect unique Q2 and t ranges that *appear* in the binning scheme for a specific xB range
static void uniqueQT_for_xB(
    const std::vector<Binning>& scheme,
    const std::pair<double,double>& xBrange,
    std::vector<std::pair<double,double>>& Q2_list,
    std::vector<std::pair<double,double>>& t_list
) {
    std::set<std::pair<double,double>> qs, ts;
    for (const auto& b : scheme) {
        if (std::make_pair(b.xBmin,b.xBmax) == xBrange) {
            qs.emplace(b.Q2min,b.Q2max);
            ts.emplace(b.tmin,b.tmax);
        }
    }
    Q2_list.assign(qs.begin(), qs.end());
    t_list.assign(ts.begin(), ts.end());
}

// Find overall indices in the global (Q2_bins, t_bins) arrays for a given range
static int findIndex(const std::pair<double,double>& range,
                     const std::vector<std::pair<double,double>>& ranges) {
    for (int i = 0; i < static_cast<int>(ranges.size()); ++i) {
        if (ranges[i] == range) return i;
    }
    return -1;
}

// Build and save canvases: one canvas per xB bin; rows=Q2 bins-in-slice, cols=t bins-in-slice.
// Each pad: contamination vs phi, with two series (+1 and -1 helicities) and error bars.
static void plotContaminationCanvases(
    const std::string& period,
    const std::map<BinKey, BinCounts>& table,
    const std::vector<Binning>& binning_scheme,
    const std::vector<std::pair<double,double>>& xB_bins,
    const std::vector<std::pair<double,double>>& Q2_bins,
    const std::vector<std::pair<double,double>>& t_bins,
    const std::string& out_dir_plots
) {
    // Prepare phi axis
    static const std::vector<double> PHI = phiCentersDeg();
    std::vector<double> X(N_PHI_BINS_PLOT), ex(N_PHI_BINS_PLOT, 0.0);
    for (int i=0;i<N_PHI_BINS_PLOT;++i) X[i] = PHI[i];

    // Ensure directory exists
    std::error_code ec;
    std::filesystem::create_directories(out_dir_plots, ec);

    // Loop over xB slices that actually appear in the scheme
    for (int ix = 0; ix < static_cast<int>(xB_bins.size()); ++ix) {
        const auto xBr = xB_bins[ix];

        // Find Q2/t grids present in this xB slice
        std::vector<std::pair<double,double>> Q2_slice, t_slice;
        uniqueQT_for_xB(binning_scheme, xBr, Q2_slice, t_slice);
        if (Q2_slice.empty() || t_slice.empty()) continue;

        const int nrows = static_cast<int>(Q2_slice.size());
        const int ncols = static_cast<int>(t_slice.size());

        // Canvas per xB slice (use two pads: a thin title pad at top and a grid pad below)
        const int w = 260*ncols + 120;
        const int h = 220*nrows + 120;

        std::ostringstream cname;
        cname << "c_contam_" << period << "_xB" << ix;
        TCanvas* c = new TCanvas(cname.str().c_str(), cname.str().c_str(), w, h);

        // Create a top title pad to avoid clipping/overlap with subpads
        TPad* pTop  = new TPad("pTop","pTop", 0.0, 0.94, 1.0, 1.0);
        pTop->SetFillStyle(0);
        pTop->SetBorderSize(0);
        pTop->Draw();

        // Main grid pad contains the divided subpads
        TPad* pGrid = new TPad("pGrid","pGrid", 0.0, 0.0, 1.0, 0.94);
        pGrid->SetFillStyle(0);
        pGrid->SetBorderSize(0);
        pGrid->Draw();
        pGrid->cd();
        pGrid->Divide(ncols, nrows, 0.0001, 0.0001);

        // Title in the dedicated top pad
        pTop->cd();
        TLatex head;
        head.SetNDC();
        head.SetTextSize(0.55);  // relative to the small top pad -> nice readable title
        head.SetTextAlign(22);   // center
        head.DrawLatex(0.5, 0.5, period.c_str());

        // Loop pads (inside pGrid)
        for (int r = 0; r < nrows; ++r) {
            int iQ2 = findIndex(Q2_slice[r], Q2_bins); if (iQ2 < 0) continue;
            for (int ccol = 0; ccol < ncols; ++ccol) {
                int itb = findIndex(t_slice[ccol], t_bins); if (itb < 0) continue;

                pGrid->cd(r*ncols + ccol + 1);
                gPad->SetGrid(1,1);
                gPad->SetTopMargin(0.08);
                gPad->SetBottomMargin(0.14);
                gPad->SetLeftMargin(0.12);
                gPad->SetRightMargin(0.06);

                // Build Y arrays (two helicities) FROM TABLE
                std::vector<double> Yp(N_PHI_BINS_PLOT, 0.0), Ym(N_PHI_BINS_PLOT, 0.0);
                std::vector<double> eYp(N_PHI_BINS_PLOT, 0.0), eYm(N_PHI_BINS_PLOT, 0.0);

                double ymax_found = 0.0;
                for (int ip = 0; ip < N_PHI_BINS_PLOT; ++ip) {
                    auto it = table.find(BinKey(ix, iQ2, itb, ip));
                    if (it == table.end()) continue;
                    const BinCounts& bc = it->second;
                    Yp[ip]  = bc.c_plus;
                    Ym[ip]  = bc.c_minus;
                    eYp[ip] = bc.c_plus_err;
                    eYm[ip] = bc.c_minus_err;
                    ymax_found = std::max(ymax_found, std::max(Yp[ip] + eYp[ip], Ym[ip] + eYm[ip]));
                }

                // y-axis range with headroom; cap at 1.0
                const double ymin = 0.0;
                const double ymax = std::min(1.0, (ymax_found > 0.0 ? ymax_found*1.25 : 0.10));

                // Graphs
                TGraphErrors* grP = new TGraphErrors(N_PHI_BINS_PLOT, X.data(), Yp.data(), ex.data(), eYp.data());
                TGraphErrors* grM = new TGraphErrors(N_PHI_BINS_PLOT, X.data(), Ym.data(), ex.data(), eYm.data());

                grP->SetTitle("");
                grM->SetTitle("");
                grP->SetMarkerStyle(24); // open circle
                grM->SetMarkerStyle(20); // filled circle
                grP->SetMarkerSize(0.9);
                grM->SetMarkerSize(0.9);
                grP->SetLineWidth(2);
                grM->SetLineWidth(2);
                grP->SetLineColor(kBlue+1);
                grP->SetMarkerColor(kBlue+1);
                grM->SetLineColor(kRed+1);
                grM->SetMarkerColor(kRed+1);

                // Draw with frame via dummy histogram axis
                TH1 *frame = gPad->DrawFrame(0.0, ymin, 360.0, ymax);
                frame->GetXaxis()->SetTitle("#phi (deg)");
                frame->GetYaxis()->SetTitle("#pi^{0} contamination");
                frame->GetXaxis()->SetNdivisions(505);
                frame->GetXaxis()->CenterTitle();
                frame->GetYaxis()->CenterTitle();

                grP->Draw("P SAME");
                grM->Draw("P SAME");

                // Fixed legend position (top-right)
                TLegend* leg = new TLegend(0.60, 0.73, 0.92, 0.92);
                leg->SetBorderSize(1);
                leg->SetLineColor(kBlack);
                leg->SetFillStyle(0);
                leg->SetTextSize(0.035);
                leg->AddEntry(grP, "helicity +1", "p");
                leg->AddEntry(grM, "helicity -1", "p");
                leg->Draw();
            }
        }

        // Save
        std::ostringstream fout;
        fout << out_dir_plots << "/plot_contamination_" << period << "_xB_" << ix << ".png";
        c->SaveAs(fout.str().c_str());

        delete c;
    }
}

// -------------------- Plot-from-JSON (public helper) --------------------
// Minimal reader for the contamination JSON written by writeContaminationJson() above.
// Loads only the fields we need for plotting.
void plot_pi0_contamination_from_json(
    const std::string& period,
    const std::vector<Binning>& binning_scheme,
    const std::string& contamination_json_path,
    const std::string& out_dir_plots
) {
    // Recover the global ranges by scanning the scheme
    auto xB_bins = uniqueRanges(binning_scheme, 'x');
    auto Q2_bins = uniqueRanges(binning_scheme, 'Q');
    auto t_bins  = uniqueRanges(binning_scheme, 't');
    if (xB_bins.empty() || Q2_bins.empty() || t_bins.empty()) {
        std::cerr << "[pi0_contam][plot-from-json][ERROR] Missing binning ranges.\n";
        return;
    }

    // Hand-rolled tiny parser (matches our write format)
    std::ifstream ifs(contamination_json_path);
    if (!ifs) {
        std::cerr << "[pi0_contam][plot-from-json][ERROR] Cannot open " << contamination_json_path << "\n";
        return;
    }
    std::string s((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    ifs.close();

    std::map<BinKey, BinCounts> table;

    size_t pos = s.find("\"bins\"");
    if (pos == std::string::npos) { std::cerr << "[plot-from-json] 'bins' not found.\n"; return; }
    pos = s.find('{', pos);
    if (pos == std::string::npos) return;
    int depth = 0; size_t i = pos;
    for (; i < s.size(); ++i) { if (s[i]=='{') depth++; else if (s[i]=='}'){ depth--; if (!depth){ ++i; break; } } }
    std::string binsObj = s.substr(pos, i-pos);

    size_t kpos = 0;
    while (true) {
        size_t keyS = binsObj.find('"', kpos);
        if (keyS == std::string::npos) break;
        size_t keyE = binsObj.find('"', keyS+1);
        if (keyE == std::string::npos) break;
        std::string key = binsObj.substr(keyS+1, keyE-keyS-1); // "(ix,iQ2,it,ip)"
        int ix, iQ2, itb, ip;
        if (std::sscanf(key.c_str(), "(%d,%d,%d,%d)", &ix, &iQ2, &itb, &ip) != 4) { kpos = keyE+1; continue; }

        size_t objS = binsObj.find('{', keyE);
        if (objS == std::string::npos) break;
        int d2=0; size_t j=objS;
        for (; j < binsObj.size(); ++j) { if (binsObj[j]=='{') d2++; else if (binsObj[j]=='}'){ d2--; if(!d2){ ++j; break; } } }
        std::string obj = binsObj.substr(objS, j-objS);

        auto findVal = [&](const std::string& path)->long long {
            size_t p = obj.find(path);
            if (p == std::string::npos) return 0;
            p = obj.find(':', p); if (p == std::string::npos) return 0;
            size_t a = p+1;
            while (a<obj.size() && std::isspace(static_cast<unsigned char>(obj[a]))) ++a;
            size_t b=a;
            while (b<obj.size() && (std::isdigit(static_cast<unsigned char>(obj[b])) || obj[b]=='-')) ++b;
            try { return std::stoll(obj.substr(a,b-a)); } catch (...) { return 0; }
        };
        auto findDouble = [&](const std::string& path)->double {
            size_t p = obj.find(path);
            if (p == std::string::npos) return 0.0;
            p = obj.find(':', p); if (p == std::string::npos) return 0.0;
            size_t a = p+1;
            while (a<obj.size() && std::isspace(static_cast<unsigned char>(obj[a]))) ++a;
            size_t b=a;
            while (b<obj.size() && (std::isdigit(static_cast<unsigned char>(obj[b])) || obj[b]=='-' || obj[b]=='+' || obj[b]=='.' || obj[b]=='e' || obj[b]=='E')) ++b;
            try { return std::stod(obj.substr(a,b-a)); } catch (...) { return 0.0; }
        };

        BinCounts bc;
        bc.N_data.plus   = findVal("\"N_data\":{\"helicity\":{\"+1\"");
        bc.N_data.minus  = findVal("\"N_data\":{\"helicity\":{\"-1\"");
        bc.N_pi0_exp.plus  = findVal("\"N_pi0_exp\":{\"helicity\":{\"+1\"");
        bc.N_pi0_exp.minus = findVal("\"N_pi0_exp\":{\"helicity\":{\"-1\"");
        bc.N_pi0_mc     = static_cast<long long>(findDouble("\"N_pi0_mc\""));
        bc.N_pi0_reco   = static_cast<long long>(findDouble("\"N_pi0_reco\""));
        bc.c_plus       = findDouble("\"contamination\":{\"+1\":{\"value\"");
        bc.c_plus_err   = findDouble("\"contamination\":{\"+1\":{\"err\"");
        bc.c_minus      = findDouble("\"contamination\":{\"-1\":{\"value\"");
        bc.c_minus_err  = findDouble("\"contamination\":{\"-1\":{\"err\"");

        table[BinKey(ix,iQ2,itb,ip)] = bc;

        kpos = j;
    }

    plotContaminationCanvases(period, table, binning_scheme, xB_bins, Q2_bins, t_bins, out_dir_plots);
}

// -------------- Combined JSON writer ----------------
static void writeCombinedContaminationJson(
    const std::string& out_path_combined,
    const std::map<std::string, std::map<BinKey, BinCounts>>& byPeriod,
    int nPhi,
    const std::vector<std::pair<double,double>>& xB_bins,
    const std::vector<std::pair<double,double>>& Q2_bins,
    const std::vector<std::pair<double,double>>& t_bins)
{
    std::ofstream ofs(out_path_combined);
    if (!ofs) { std::cerr << "[pi0_contam][ERROR] Cannot open combined output " << out_path_combined << "\n"; return; }
    ofs << std::fixed << std::setprecision(8);
    ofs << "{\n";
    ofs << "  \"binning_meta\": {\"phi_bins\": " << nPhi
        << ", \"xB_bins\": " << xB_bins.size()
        << ", \"Q2_bins\": " << Q2_bins.size()
        << ", \"t_bins\": "  << t_bins.size()  << "},\n";
    ofs << "  \"periods\": {\n";

    bool firstP = true;
    for (const auto& pkv : byPeriod) {
        if (!firstP) ofs << ",\n";
        firstP = false;
        ofs << "    \"" << pkv.first << "\": {\n";
        ofs << "      \"bins\": {\n";
        bool firstB = true;
        for (const auto& kv : pkv.second) {
            if (!firstB) ofs << ",\n";
            firstB = false;
            const auto& bc = kv.second;
            ofs << "        \"" << keyStr(kv.first) << "\": {"
                << "\"N_data\":{\"helicity\":{\"+1\":" << bc.N_data.plus
                << ",\"-1\":" << bc.N_data.minus << "},\"total\":" << (bc.N_data.plus+bc.N_data.minus) << "}"
                << ",\"N_pi0_exp\":{\"helicity\":{\"+1\":" << bc.N_pi0_exp.plus
                << ",\"-1\":" << bc.N_pi0_exp.minus << "},\"total\":" << (bc.N_pi0_exp.plus+bc.N_pi0_exp.minus) << "}"
                << ",\"N_pi0_mc\":"   << bc.N_pi0_mc
                << ",\"N_pi0_reco\":" << bc.N_pi0_reco
                << ",\"contamination\":{"
                << "\"+1\":{\"value\":" << bc.c_plus  << ",\"err\":" << bc.c_plus_err  << "},"
                << "\"-1\":{\"value\":" << bc.c_minus << ",\"err\":" << bc.c_minus_err << "}"
                << "}"
                << "}";
        }
        ofs << "\n      }\n"; // bins
        ofs << "    }";
    }

    ofs << "\n  }\n}\n";
    ofs.close();
    std::cout << "[pi0_contam] Wrote combined " << out_path_combined << std::endl;
}

// ---------------- core ----------------
void compute_pi0_contamination_helicity(
    const std::vector<std::string>& periods,
    const std::vector<std::string>& topologies,
    const std::vector<Binning>& binning_scheme,
    const std::map<std::string, TTree*>& dvcsDataTrees,
    const std::map<std::string, TTree*>& eppi0DataTrees,
    const std::map<std::string, TTree*>& eppi0RecMcTrees,
    const std::map<std::string, TTree*>& eppi0BkgTrees,
    const std::string& combined_cuts_json,
    const std::string& out_root_dir
) {
    namespace fs = std::filesystem;

    // Directories (per your conventions):
    const fs::path root(out_root_dir); // expected "output"
    const fs::path jsons_dir = root / "jsons" / "contamination";             // per-period JSONs
    const fs::path combined_json_path = root / "jsons" / "pi0_contamination_combined.json"; // combined JSON
    const fs::path plots_root = root / "contamination_plots";                 // already created by make_dirs

    std::error_code ec;
    fs::create_directories(jsons_dir, ec);

    const auto xB_bins = uniqueRanges(binning_scheme, 'x');
    const auto Q2_bins = uniqueRanges(binning_scheme, 'Q');
    const auto t_bins  = uniqueRanges(binning_scheme, 't');
    if (xB_bins.empty() || Q2_bins.empty() || t_bins.empty()) {
        std::cerr << "[pi0_contam][ERROR] Missing binning ranges." << std::endl;
        return;
    }

    // Load exclusivity 3σ cuts (data) per DVCS period/topology
    PeriodTopoCuts cuts;
    if (!loadCombinedCuts(combined_cuts_json, cuts)) {
        std::cerr << "[pi0_contam][WARN] No combined cuts loaded; proceeding without 3σ cuts." << std::endl;
    }

    auto dvcsCutsKey = [&](const std::string& runTag, const std::string& topoKey)->std::string {
        std::string cap = runTag;
        if (!cap.empty()) cap[0] = std::toupper(cap[0]);
        for (size_t i=0;i+1<cap.size();++i) if (cap[i]=='_' && i+1<cap.size()) cap[i+1]=std::toupper(cap[i+1]);
        return std::string("DVCS_") + cap + "_" + topoKey;
    };

    // Keep a copy of all per-period tables to write the combined file
    std::map<std::string, std::map<BinKey, BinCounts>> allPeriods;

    // For each DVCS period, build counts
    for (const auto& period : periods) {
        const std::string runTag = periodToRunTagKey(period);       // e.g. fa18_inb
        const std::string key_dvcs = runTag;                        // dvcs data key
        const std::string key_pi0_data = runTag + "_eppi0";         // eppi0 data
        const std::string key_pi0_reco = runTag + "_rec_mc";        // eppi0 reco MC
        const std::string key_pi0_bkg  = runTag + "_bkg";           // eppi0 mis-ID MC (to DVCS)

        auto itDVCS = dvcsDataTrees.find(key_dvcs);
        if (itDVCS == dvcsDataTrees.end() || !itDVCS->second) {
            std::cerr << "[pi0_contam][WARN] No DVCS data tree for '" << period << "' (key '" << key_dvcs << "'). Skipping." << std::endl;
            continue;
        }
        TTree* t_dvcs = itDVCS->second;

        TTree* t_pi0_data = nullptr;
        auto itPi0Data = eppi0DataTrees.find(key_pi0_data);
        if (itPi0Data != eppi0DataTrees.end()) t_pi0_data = itPi0Data->second;

        TTree* t_pi0_reco = nullptr;
        auto itPi0Reco = eppi0RecMcTrees.find(key_pi0_reco);
        if (itPi0Reco != eppi0RecMcTrees.end()) t_pi0_reco = itPi0Reco->second;

        TTree* t_pi0_bkg = nullptr;
        auto itPi0Bkg = eppi0BkgTrees.find(key_pi0_bkg);
        if (itPi0Bkg != eppi0BkgTrees.end()) t_pi0_bkg = itPi0Bkg->second;

        if (!t_pi0_bkg || !t_pi0_reco) {
            std::cerr << "[pi0_contam][WARN] Missing π0 MC for '" << period << "'; bkg=" << (t_pi0_bkg?"ok":"none")
                      << " reco=" << (t_pi0_reco?"ok":"none") << ". Continuing with what is available." << std::endl;
        }
        if (!t_pi0_data) {
            std::cerr << "[pi0_contam][WARN] Missing eppi0 DATA for '" << period << "'. Continuing with zeros for N_pi0_exp." << std::endl;
        }

        // Accumulator per bin
        std::map<BinKey, BinCounts> counts;

        // ---- DVCS data (helicity-resolved) ----
        {
            BranchBinderDVCS b; b.bind(t_dvcs);
            if (!b.readyCuts() || !b.readyBins()) {
                std::cerr << "[pi0_contam][WARN] DVCS tree missing branches for '" << period << "'. Skipping DVCS loop." << std::endl;
            } else {
                const Long64_t nent = t_dvcs->GetEntries();
                for (Long64_t i=0;i<nent;++i) {
                    t_dvcs->GetEntry(i);
                    if (!applyKinematicCuts_simple(b.t1,b.open_angle_ep2,b.pTmiss)) continue;
                    if (b.helicity!=+1 && b.helicity!=-1) continue;

                    // choose topology for this event
                    std::string usedTopoKey;
                    for (const auto& topoStr : topologies) {
                        if (passesTopology_simple(b.detector1,b.detector2,topoStr)) { usedTopoKey = topoToKey(topoStr); break; }
                    }
                    if (usedTopoKey.empty()) continue;

                    // 3σ cuts (if available)
                    VarCutMap topoCuts;
                    auto itCut = cuts.find(dvcsCutsKey(runTag, usedTopoKey));
                    if (itCut != cuts.end()) topoCuts = itCut->second;
                    if (!topoCuts.empty()) {
                        if (!passes3SigmaCuts(topoCuts, b.cutVals())) continue;
                    }

                    // binning
                    double xB=b.x, Q2=b.Q2, tt=std::fabs(b.t1), phi=b.phi();
                    if (!std::isfinite(xB)||!std::isfinite(Q2)||!std::isfinite(tt)||!std::isfinite(phi)) continue;
                    int ix=findBin(xB,xB_bins), iQ=findBin(Q2,Q2_bins), itb=findBin(tt,t_bins), ip=phiToBin(phi);
                    if (ix<0||iQ<0||itb<0||ip<0||ip>=N_PHI_BINS) continue;

                    BinKey key(ix,iQ,itb,ip);
                    if (b.helicity==+1) counts[key].N_data.plus++; else counts[key].N_data.minus++;
                }
            }
        }

        // ---- π0 background MC mis-ID to DVCS (no helicity) -> N_pi0_mc ----
        if (t_pi0_bkg) {
            BranchBinderEPPI0MC b; b.bind(t_pi0_bkg);
            if (!b.readyCuts() || !b.readyBins()) {
                std::cerr << "[pi0_contam][WARN] eppi0_bkg tree missing branches for '" << period << "'. Skipping." << std::endl;
            } else {
                const Long64_t nent = t_pi0_bkg->GetEntries();
                for (Long64_t i=0;i<nent;++i) {
                    t_pi0_bkg->GetEntry(i);
                    if (!applyKinematicCuts_simple(b.t1,b.open_angle_ep2,b.pTmiss)) continue;

                    bool match=false; std::string usedTopoKey;
                    for (const auto& topoStr : topologies) {
                        if (passesTopology_simple(b.detector1,b.detector2,topoStr)) { usedTopoKey = topoToKey(topoStr); match=true; break; }
                    }
                    if (!match) continue;

                    // 3σ cuts under DVCS hypothesis
                    VarCutMap topoCuts;
                    auto itCut = cuts.find(dvcsCutsKey(runTag, usedTopoKey));
                    if (itCut != cuts.end()) topoCuts = itCut->second;
                    if (!topoCuts.empty()) {
                        if (!passes3SigmaCuts(topoCuts, b.cutValsForDVCS())) continue;
                    }

                    // binning
                    double xB=b.x, Q2=b.Q2, tt=std::fabs(b.t1), phi=b.phi();
                    if (!std::isfinite(xB)||!std::isfinite(Q2)||!std::isfinite(tt)||!std::isfinite(phi)) continue;
                    int ix=findBin(xB,xB_bins), iQ=findBin(Q2,Q2_bins), itb=findBin(tt,t_bins), ip=phiToBin(phi);
                    if (ix<0||iQ<0||itb<0||ip<0||ip>=N_PHI_BINS) continue;
                    counts[BinKey(ix,iQ,itb,ip)].N_pi0_mc++;
                }
            }
        }

        // ---- π0 reco MC (no helicity) -> N_pi0_reco ----
        if (t_pi0_reco) {
            BranchBinderEPPI0MC b; b.bind(t_pi0_reco);
            if (!b.readyCuts() || !b.readyBins()) {
                std::cerr << "[pi0_contam][WARN] eppi0 reco MC tree missing branches for '" << period << "'. Skipping." << std::endl;
            } else {
                const Long64_t nent = t_pi0_reco->GetEntries();
                for (Long64_t i=0;i<nent;++i) {
                    t_pi0_reco->GetEntry(i);
                    if (!applyKinematicCuts_simple(b.t1,b.open_angle_ep2,b.pTmiss)) continue;

                    bool match=false;
                    for (const auto& topoStr : topologies) {
                        if (passesTopology_simple(b.detector1,b.detector2,topoStr)) { match=true; break; }
                    }
                    if (!match) continue;

                    // binning
                    double xB=b.x, Q2=b.Q2, tt=std::fabs(b.t1), phi=b.phi();
                    if (!std::isfinite(xB)||!std::isfinite(Q2)||!std::isfinite(tt)||!std::isfinite(phi)) continue;
                    int ix=findBin(xB,xB_bins), iQ=findBin(Q2,Q2_bins), itb=findBin(tt,t_bins), ip=phiToBin(phi);
                    if (ix<0||iQ<0||itb<0||ip<0||ip>=N_PHI_BINS) continue;
                    counts[BinKey(ix,iQ,itb,ip)].N_pi0_reco++;
                }
            }
        }

        // ---- eppi0 experimental data (helicity-resolved) -> N_pi0_exp^± ----
        if (t_pi0_data) {
            BranchBinderEPPI0Data b; b.bind(t_pi0_data);
            if (!b.readyCuts() || !b.readyBins()) {
                std::cerr << "[pi0_contam][WARN] eppi0 DATA tree missing branches for '" << period << "'. Skipping eppi0 DATA." << std::endl;
            } else {
                const Long64_t nent = t_pi0_data->GetEntries();
                for (Long64_t i=0;i<nent;++i) {
                    t_pi0_data->GetEntry(i);
                    if (!applyKinematicCuts_simple(b.t1,b.open_angle_ep2,b.pTmiss)) continue;
                    if (b.helicity!=+1 && b.helicity!=-1) continue;

                    bool match=false;
                    for (const auto& topoStr : topologies) {
                        if (passesTopology_simple(b.detector1,b.detector2,topoStr)) { match=true; break; }
                    }
                    if (!match) continue;

                    // binning
                    double xB=b.x, Q2=b.Q2, tt=std::fabs(b.t1), phi=b.phi();
                    if (!std::isfinite(xB)||!std::isfinite(Q2)||!std::isfinite(tt)||!std::isfinite(phi)) continue;
                    int ix=findBin(xB,xB_bins), iQ=findBin(Q2,Q2_bins), itb=findBin(tt,t_bins), ip=phiToBin(phi);
                    if (ix<0||iQ<0||itb<0||ip<0||ip>=N_PHI_BINS) continue;

                    BinKey key(ix,iQ,itb,ip);
                    if (b.helicity==+1) counts[key].N_pi0_exp.plus++; else counts[key].N_pi0_exp.minus++;
                }
            }
        }

        // ---- compute contamination per helicity ----
        for (auto& kv : counts) {
            BinCounts& bc = kv.second;

            auto compute_one = [&](long long N_pi0_exp_h, long long N_data_h, double& c, double& c_err){
                if (N_data_h<=0 || bc.N_pi0_reco<=0 || bc.N_pi0_mc<=0) { c=0.0; c_err=0.0; return; }
                double ratio = static_cast<double>(N_pi0_exp_h) / static_cast<double>(bc.N_pi0_reco);
                double cval  = static_cast<double>(bc.N_pi0_mc) * ratio / static_cast<double>(N_data_h);

                auto rel = [](long long n)->double { return (n>0)? 1.0/std::sqrt(static_cast<double>(n)) : 0.0; };
                double rel_mc   = rel(bc.N_pi0_mc);
                double rel_exp  = rel(N_pi0_exp_h);
                double rel_reco = rel(bc.N_pi0_reco);
                double rel_data = rel(N_data_h);

                double rel_ratio = std::sqrt(rel_exp*rel_exp + rel_reco*rel_reco);
                double rel_tot   = std::sqrt(rel_mc*rel_mc + rel_ratio*rel_ratio + rel_data*rel_data);
                c = cval;
                c_err = cval * rel_tot;
            };

            compute_one(bc.N_pi0_exp.plus,  bc.N_data.plus,  bc.c_plus,  bc.c_plus_err);
            compute_one(bc.N_pi0_exp.minus, bc.N_data.minus, bc.c_minus, bc.c_minus_err);
        }

        // ---- write JSON for this period ----
        const std::string out_file = (jsons_dir / ("contamination_" + period + ".json")).string();
        writeContaminationJson(out_file, counts, N_PHI_BINS, xB_bins, Q2_bins, t_bins);

        // Plot canvases for this period (to output/contamination_plots/<runTag>/...)
        if (ENABLE_PI0_CONTAMINATION_PLOTS) {
            const fs::path period_plot_dir = plots_root / runTag;
            std::error_code ec_plot;
            fs::create_directories(period_plot_dir, ec_plot);
            plotContaminationCanvases(period, counts, binning_scheme, xB_bins, Q2_bins, t_bins, period_plot_dir.string());
        }

        // ---- optional copy for Fa18_inb_supp (JSON only; no plots) ----
        if (COPY_CONTAM_TO_FA18_INB_SUPP && runTag == "fa18_inb") {
            const std::string supp_period = "DVCS_Fa18_inb_supp";
            const std::string out_copy = (jsons_dir / ("contamination_" + supp_period + ".json")).string();
            std::error_code ec_copy;
            fs::copy_file(out_file, out_copy, fs::copy_options::overwrite_existing, ec_copy);
            if (ec_copy) std::cerr << "[pi0_contam][WARN] Could not copy to Fa18_inb_supp JSON: " << ec_copy.message() << std::endl;
            else         std::cout << "[pi0_contam] Also wrote (copy) " << out_copy << std::endl;

            // Do NOT plot for Fa18_inb_supp to avoid duplicate identical plots
            allPeriods[supp_period] = counts; // include in combined
        }

        // Keep for combined output
        allPeriods[period] = counts;
    } // periods

    // ---- write combined JSON with everything ----
    writeCombinedContaminationJson(combined_json_path.string(), allPeriods, N_PHI_BINS, xB_bins, Q2_bins, t_bins);
}

--

Here is total_counts.cpp:

#include "total_counts.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

// ---------------- simple helpers ----------------
static inline std::string toLower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(),
        [](unsigned char c){ return static_cast<char>(std::tolower(c)); });
    return s;
}

static inline std::string trim(const std::string& s) {
    size_t a = s.find_first_not_of(" \t\r\n");
    size_t b = s.find_last_not_of(" \t\r\n");
    if (a == std::string::npos) return "";
    return s.substr(a, b - a + 1);
}

// Period like "DVCS_Fa18_inb" -> key "fa18_inb"; "DVCS_Sp18_out" -> "sp18_out".
static std::string periodToRunTagKey(const std::string& period) {
    auto pos = period.find('_');
    if (pos == std::string::npos || pos + 1 >= period.size()) return toLower(period);
    return toLower(period.substr(pos + 1));
}

// Map run-tag "fa18_inb" -> "DVCS_Fa18_inb"
static inline std::string dvcsPeriodName(const std::string& runTag) {
    std::string cap = runTag;
    if (!cap.empty()) cap[0] = static_cast<char>(std::toupper(static_cast<unsigned char>(cap[0])));
    for (size_t i = 0; i + 1 < cap.size(); ++i) {
        if (cap[i] == '_' && i + 1 < cap.size()) {
            cap[i + 1] = static_cast<char>(std::toupper(static_cast<unsigned char>(cap[i + 1])));
        }
    }
    return std::string("DVCS_") + cap;
}

static constexpr int    N_PHI_BINS = 12;
static constexpr double TWO_PI     = 2.0 * M_PI;

// wrap phi to [0, 2pi)
static inline double wrapToTwoPi(double phi) {
    if (!std::isfinite(phi)) return 0.0;
    double w = std::fmod(phi, TWO_PI);
    if (w < 0.0) w += TWO_PI;
    if (w >= TWO_PI) w = std::nextafter(TWO_PI, 0.0);
    return w;
}

static inline int phiToBin(double phi) {
    double w = wrapToTwoPi(phi);
    double width = TWO_PI / static_cast<double>(N_PHI_BINS);
    int idx = static_cast<int>(std::floor(w / width));
    if (idx < 0) idx = 0;
    if (idx >= N_PHI_BINS) idx = N_PHI_BINS - 1;
    return idx;
}

// Build unique (min,max) bins from the scheme for a coordinate (x, Q, t).
static std::vector<std::pair<double,double>> uniqueRanges(const std::vector<Binning>& scheme, char which) {
    std::set<std::pair<double,double>> s;
    for (const auto& b : scheme) {
        if (which == 'x') s.emplace(b.xBmin, b.xBmax);
        else if (which == 'Q') s.emplace(b.Q2min, b.Q2max);
        else if (which == 't') s.emplace(b.tmin, b.tmax);
    }
    return std::vector<std::pair<double,double>>(s.begin(), s.end());
}

static int findBin(double v, const std::vector<std::pair<double,double>>& ranges) {
    for (int i = 0; i < static_cast<int>(ranges.size()); ++i)
        if (v >= ranges[i].first && v < ranges[i].second) return i;
    return -1;
}

// ------------- topology + universal kinematic cuts -------------
static inline bool passesTopology_simple(int detector1, int detector2, const std::string& topoStr) {
    if (topoStr == "(FD,FD)") return (detector1 == 1 && detector2 == 1);
    if (topoStr == "(CD,FD)") return (detector1 == 2 && detector2 == 1);
    if (topoStr == "(CD,FT)") return (detector1 == 2 && detector2 == 0);
    return false;
}

static inline std::string topoToKey(const std::string& topoStr) {
    if (topoStr == "(FD,FD)") return "FD_FD";
    if (topoStr == "(CD,FD)") return "CD_FD";
    if (topoStr == "(CD,FT)") return "CD_FT";
    return "FD_FD";
}

static inline bool applyKinematicCuts_simple(double t1, double open_angle_ep2, double pTmiss) {
    if (open_angle_ep2 <= 5.0) return false;
    if ((-t1) > 1.0)          return false;
    if (pTmiss > 0.20)        return false;
    return true;
}

// ------------- exclusivity 3sigma cuts: data-only stats -------------
struct Stats { double mean = 0.0; double std = 0.0; };
using VarCutMap = std::map<std::string, Stats>;       // var -> {mean,std}
using PeriodTopoCuts = std::map<std::string, VarCutMap>; // "DVCS_Sp18_inb_FD_FD" -> data cuts

// We expect JSON written by exclusivity_cuts.cpp -> writeCombinedJson().
// Very small hand-rolled parser for that specific structure.
// Only extracts the "data" section for each key.
static bool loadCombinedCuts(const std::string& path, PeriodTopoCuts& out) {
    std::ifstream ifs(path);
    if (!ifs) {
        std::cerr << "[total_counts][ERROR] Cannot open cuts JSON: " << path << std::endl;
        return false;
    }
    std::string s((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    ifs.close();

    size_t pos = 0;
    while (true) {
        size_t keyStart = s.find('"', pos);
        if (keyStart == std::string::npos) break;
        size_t keyEnd = s.find('"', keyStart + 1);
        if (keyEnd == std::string::npos) break;
        std::string key = s.substr(keyStart + 1, keyEnd - keyStart - 1); // e.g., DVCS_Sp18_inb_FD_FD

        size_t dataPos = s.find("\"data\"", keyEnd);
        if (dataPos == std::string::npos) { pos = keyEnd + 1; continue; }
        size_t braceStart = s.find('{', dataPos);
        if (braceStart == std::string::npos) { pos = keyEnd + 1; continue; }

        int depth = 0; size_t i = braceStart;
        for (; i < s.size(); ++i) {
            if (s[i] == '{') depth++;
            else if (s[i] == '}') {
                depth--;
                if (depth == 0) { ++i; break; }
            }
        }
        if (depth != 0) { pos = keyEnd + 1; continue; }
        std::string dataObj = s.substr(braceStart, i - braceStart);

        VarCutMap cuts;
        size_t vpos = 0;
        while (true) {
            size_t vKeyS = dataObj.find('"', vpos);
            if (vKeyS == std::string::npos) break;
            size_t vKeyE = dataObj.find('"', vKeyS + 1);
            if (vKeyE == std::string::npos) break;
            std::string var = dataObj.substr(vKeyS + 1, vKeyE - vKeyS - 1);

            size_t meanPos = dataObj.find("\"mean\"", vKeyE);
            size_t stdPos  = dataObj.find("\"std\"",  vKeyE);
            if (meanPos == std::string::npos || stdPos == std::string::npos) { vpos = vKeyE + 1; continue; }

            auto readNumberAfterColon = [&](size_t from)->double {
                size_t colon = dataObj.find(':', from);
                if (colon == std::string::npos) return 0.0;
                size_t j = colon + 1;
                while (j < dataObj.size() && std::isspace(static_cast<unsigned char>(dataObj[j]))) ++j;
                size_t k = j;
                while (k < dataObj.size() &&
                       (std::isdigit(static_cast<unsigned char>(dataObj[k])) ||
                        dataObj[k] == '-' || dataObj[k] == '+' || dataObj[k] == '.' ||
                        dataObj[k] == 'e' || dataObj[k] == 'E')) ++k;
                std::string num = dataObj.substr(j, k - j);
                try { return std::stod(num); } catch (...) { return 0.0; }
            };

            double m = readNumberAfterColon(meanPos);
            double sdev = readNumberAfterColon(stdPos);
            cuts[var] = Stats{m, sdev};

            vpos = vKeyE + 1;
        }

        if (key.rfind("DVCS_", 0) == 0) out[key] = cuts;

        pos = keyEnd + 1;
    }

    if (out.empty()) {
        std::cerr << "[total_counts][WARN] No DVCS cuts found in " << path << std::endl;
        return false;
    }
    return true;
}

static inline bool within3Sigma(double val, const Stats& s) {
    return (val >= s.mean - 3.0 * s.std) && (val <= s.mean + 3.0 * s.std);
}

static bool passes3SigmaCuts(const VarCutMap& cuts, const std::map<std::string,double>& values) {
    for (const auto& kv : cuts) {
        auto it = values.find(kv.first);
        if (it == values.end()) continue;
        if (!within3Sigma(it->second, kv.second)) return false;
    }
    return true;
}

// ------------- Branch binder (includes helicity and bin vars) -------------
struct BranchBinderTC {
    // topology helpers
    int detector1 = 0; bool has_detector1 = false;
    int detector2 = 0; bool has_detector2 = false;

    // kinematic cuts
    double t1 = 0.0;             bool has_t1 = false;
    double open_angle_ep2 = 0.0; bool has_open_angle_ep2 = false;
    double pTmiss = 0.0;         bool has_pTmiss = false;

    // binning variables
    double x = 0.0;     bool has_x = false;   // Bjorken x
    double Q2 = 0.0;    bool has_Q2 = false;
    double phi2 = 0.0;  bool has_phi2 = false;
    double Delta_phi = 0.0; bool has_Delta_phi = false;

    // exclusivity variables (for 3sigma checks)
    double Emiss2 = 0.0;            bool has_Emiss2 = false;
    double Mx2 = 0.0;               bool has_Mx2 = false;
    double Mx2_1 = 0.0;             bool has_Mx2_1 = false;
    double Mx2_2 = 0.0;             bool has_Mx2_2 = false;
    double theta_gamma_gamma = 0.0; bool has_theta_gamma_gamma = false;
    double xF = 0.0;                bool has_xF = false;

    int helicity = 0; bool has_helicity = false;

    void bind(TTree* t) {
        if (!t) return;
        auto bindI = [&](const char* n, int* a, bool& f){ if (t->GetBranch(n)) { t->SetBranchAddress(n, a); f = true; } };
        auto bindD = [&](const char* n, double* a, bool& f){ if (t->GetBranch(n)) { t->SetBranchAddress(n, a); f = true; } };

        bindI("detector1", &detector1, has_detector1);
        bindI("detector2", &detector2, has_detector2);
        bindI("helicity",  &helicity,  has_helicity);

        bindD("t1", &t1, has_t1);
        bindD("open_angle_ep2", &open_angle_ep2, has_open_angle_ep2);
        bindD("pTmiss", &pTmiss, has_pTmiss);

        bindD("x", &x, has_x);
        bindD("Q2", &Q2, has_Q2);
        bindD("phi2", &phi2, has_phi2);
        bindD("Delta_phi", &Delta_phi, has_Delta_phi);

        bindD("Emiss2", &Emiss2, has_Emiss2);
        bindD("Mx2", &Mx2, has_Mx2);
        bindD("Mx2_1", &Mx2_1, has_Mx2_1);
        bindD("Mx2_2", &Mx2_2, has_Mx2_2);
        bindD("theta_gamma_gamma", &theta_gamma_gamma, has_theta_gamma_gamma);
        bindD("xF", &xF, has_xF);
    }

    double phi() const {
        return has_phi2 ? phi2 :
               (has_Delta_phi ? Delta_phi : std::numeric_limits<double>::quiet_NaN());
    }

    bool readyForCuts() const {
        return has_detector1 && has_detector2 && has_t1 && has_open_angle_ep2 && has_pTmiss && has_helicity;
    }

    bool readyForBinning() const {
        return has_x && has_Q2 && (has_phi2 || has_Delta_phi);
    }

    std::map<std::string,double> valueMapForCuts() const {
        std::map<std::string,double> m;
        if (has_Delta_phi) m["Delta_phi"] = Delta_phi;
        if (has_theta_gamma_gamma) m["theta_gamma_gamma"] = theta_gamma_gamma;
        if (has_pTmiss) m["pTmiss"] = pTmiss;
        if (has_xF) m["xF"] = xF;
        if (has_Emiss2) m["Emiss2"] = Emiss2;
        if (has_Mx2) m["Mx2"] = Mx2;
        if (has_Mx2_1) m["Mx2_1"] = Mx2_1;
        if (has_Mx2_2) m["Mx2_2"] = Mx2_2;
        return m;
    }
};

// ------------- counting core -------------
struct HelCounts { long long plus = 0; long long minus = 0; };
using BinKey = std::tuple<int,int,int,int>; // (ix, iQ2, it, ip)

// group -> binKey -> HelCounts
using GroupCounts = std::map<std::string, std::map<BinKey, HelCounts>>;

static inline std::string tupleKeyStr(const BinKey& k) {
    int ix, iQ2, it, ip; std::tie(ix, iQ2, it, ip) = k;
    std::ostringstream os; os << "(" << ix << "," << iQ2 << "," << it << "," << ip << ")";
    return os.str();
}

// Individual per-period group name should be DVCS_* to match downstream readers
static inline std::string indivGroupName(const std::string& runTag) {
    return dvcsPeriodName(runTag); // e.g., "DVCS_Fa18_inb"
}

// Existing combined group names (kept for backward-compat)
static inline std::string combFa18() { return "Fa18_inb_out_combined"; }
static inline std::string combFour() { return "Sp18_and_Fa18_inb_out_combined"; }

// New requested combined group names
static inline std::string grpFa18Sum() { return "fa18_sum"; }
static inline std::string grpTenSixSum() { return "10.6_sum"; }

// Membership checks for combined groups
static inline bool inFa18(const std::string& runTag) {
    return (runTag == "fa18_inb" || runTag == "fa18_out");
}
static inline bool inFour(const std::string& runTag) {
    return (runTag == "sp18_inb" || runTag == "sp18_out" || runTag == "fa18_inb" || runTag == "fa18_out");
}
static inline bool inFa18Sum(const std::string& runTag) {
    return (runTag == "fa18_inb" || runTag == "fa18_out");
}
static inline bool inTenSixSum(const std::string& runTag) {
    return (runTag == "sp18_inb" || runTag == "sp18_out" || runTag == "fa18_inb" || runTag == "fa18_inb_supp" || runTag == "fa18_out");
}

// JSON writer
static void writeCountsJson(const std::string& outPath,
                            const GroupCounts& groups,
                            int nPhi,
                            const std::vector<std::pair<double,double>>& xB_bins,
                            const std::vector<std::pair<double,double>>& Q2_bins,
                            const std::vector<std::pair<double,double>>& t_bins)
{
    std::ofstream ofs(outPath);
    if (!ofs) {
        std::cerr << "[total_counts][ERROR] Cannot open output JSON: " << outPath << std::endl;
        return;
    }
    ofs << std::fixed << std::setprecision(8);

    ofs << "{\n";
    ofs << "  \"binning_meta\": {\"phi_bins\": " << nPhi
        << ", \"xB_bins\": " << xB_bins.size()
        << ", \"Q2_bins\": " << Q2_bins.size()
        << ", \"t_bins\": "  << t_bins.size()  << "},\n";
    ofs << "  \"groups\": {\n";

    bool firstGroup = true;
    for (const auto& gkv : groups) {
        if (!firstGroup) ofs << ",\n";
        firstGroup = false;

        ofs << "    \"" << gkv.first << "\": {\n";
        ofs << "      \"bins\": {\n";

        bool firstBin = true;
        for (const auto& bkv : gkv.second) {
            if (!firstBin) ofs << ",\n";
            firstBin = false;

            std::string keyStr = tupleKeyStr(bkv.first);
            const HelCounts& hc = bkv.second;
            long long total = hc.plus + hc.minus;

            ofs << "        \"" << keyStr << "\": {\"helicity\": {\"+1\": " << hc.plus
                << ", \"-1\": " << hc.minus << "}, \"total\": " << total << "}";
        }

        ofs << "\n      }\n";   // end bins
        ofs << "    }";         // end group
    }

    ofs << "\n  }\n}\n";
    ofs.close();
    std::cout << "[total_counts] Wrote " << outPath << std::endl;
}

void compute_total_counts(const std::vector<std::string>& periods,
                          const std::vector<std::string>& topologies,
                          const std::vector<Binning>& binning_scheme,
                          const std::map<std::string, TTree*>& dataTrees,
                          const std::string& cuts_json_path,
                          const std::string& output_json_path)
{
    // Build bin edges
    const auto xB_bins = uniqueRanges(binning_scheme, 'x');
    const auto Q2_bins = uniqueRanges(binning_scheme, 'Q');
    const auto t_bins  = uniqueRanges(binning_scheme, 't');
    if (xB_bins.empty() || Q2_bins.empty() || t_bins.empty()) {
        std::cerr << "[total_counts][ERROR] Missing xB/Q2/t bins. Aborting." << std::endl;
        return;
    }

    // Load exclusivity 3sigma cuts (data) from combined JSON
    PeriodTopoCuts cuts;
    bool ok = loadCombinedCuts(cuts_json_path, cuts);
    if (!ok) std::cerr << "[total_counts][WARN] Proceeding, but no cuts were loaded." << std::endl;

    // Set up output container
    GroupCounts outCounts;

    // Loop over each requested period
    for (const std::string& period : periods) {
        const std::string runTag = periodToRunTagKey(period); // e.g., "sp18_inb"
        auto itTree = dataTrees.find(runTag);
        if (itTree == dataTrees.end() || !itTree->second) {
            std::cerr << "[total_counts][WARN] No data tree for '" << period
                      << "' (key '" << runTag << "'). Skipping." << std::endl;
            continue;
        }
        TTree* t = itTree->second;

        // Bind branches once per tree
        BranchBinderTC b;
        b.bind(t);
        if (!b.readyForCuts() || !b.readyForBinning()) {
            std::cerr << "[total_counts][WARN] Missing required branches for '" << period << "'. Skipping." << std::endl;
            continue;
        }

        const Long64_t nent = t->GetEntries();

        // Build DVCS key prefix used in exclusivity combined JSON
        std::string cap = runTag;
        if (!cap.empty()) cap[0] = static_cast<char>(std::toupper(static_cast<unsigned char>(cap[0])));
        for (size_t i = 0; i + 1 < cap.size(); ++i) {
            if (cap[i] == '_' && i + 1 < cap.size()) {
                cap[i + 1] = static_cast<char>(std::toupper(static_cast<unsigned char>(cap[i + 1])));
            }
        }
        const std::string periodCode = std::string("DVCS_") + cap;

        // Event loop
        for (Long64_t i = 0; i < nent; ++i) {
            t->GetEntry(i);

            if (!applyKinematicCuts_simple(b.t1, b.open_angle_ep2, b.pTmiss)) continue;
            if (b.helicity != +1 && b.helicity != -1) continue;

            // Determine topology for this event
            std::string usedTopoKey;
            for (const auto& topoStr : topologies) {
                if (passesTopology_simple(b.detector1, b.detector2, topoStr)) {
                    usedTopoKey = topoToKey(topoStr);
                    break;
                }
            }
            if (usedTopoKey.empty()) continue;

            // Load cuts for this period+topology if available
            VarCutMap topoCuts;
            auto itCut = cuts.find(periodCode + "_" + usedTopoKey);
            if (itCut != cuts.end()) topoCuts = itCut->second;

            // Apply 3sigma cuts
            if (!topoCuts.empty()) {
                auto vals = b.valueMapForCuts();
                if (!passes3SigmaCuts(topoCuts, vals)) continue;
            }

            // Binning
            double xB = b.x, Q2 = b.Q2, tt = std::fabs(b.t1), phi = b.phi();
            if (!std::isfinite(xB) || !std::isfinite(Q2) || !std::isfinite(tt) || !std::isfinite(phi)) continue;

            int ix  = findBin(xB, xB_bins);
            int iQ2 = findBin(Q2, Q2_bins);
            int itb = findBin(tt,  t_bins);
            int ip  = phiToBin(phi);
            if (ix < 0 || iQ2 < 0 || itb < 0 || ip < 0 || ip >= N_PHI_BINS) continue;

            BinKey key(ix, iQ2, itb, ip);

            // Individual DVCS_* group
            HelCounts& hc_ind = outCounts[indivGroupName(runTag)][key];
            if (b.helicity == +1) hc_ind.plus++; else hc_ind.minus++;

            // Existing combined groups (kept)
            if (inFa18(runTag)) {
                HelCounts& hc_fa = outCounts[combFa18()][key];
                if (b.helicity == +1) hc_fa.plus++; else hc_fa.minus++;
            }
            if (inFour(runTag)) {
                HelCounts& hc_4 = outCounts[combFour()][key];
                if (b.helicity == +1) hc_4.plus++; else hc_4.minus++;
            }

            // New requested summed groups
            if (inFa18Sum(runTag)) {
                HelCounts& hc_fs = outCounts[grpFa18Sum()][key];
                if (b.helicity == +1) hc_fs.plus++; else hc_fs.minus++;
            }
            if (inTenSixSum(runTag)) {
                HelCounts& hc_106 = outCounts[grpTenSixSum()][key];
                if (b.helicity == +1) hc_106.plus++; else hc_106.minus++;
            }
        }

        std::cout << "[total_counts] Counted period " << period
                  << " using cuts key '" << periodCode << "_*'" << std::endl;
    }

    // Serialize JSON
    writeCountsJson(output_json_path, outCounts, N_PHI_BINS, xB_bins, Q2_bins, t_bins);
}

--

here is pi0_corrected_counts.cpp:

// pi0_corrected_counts.cpp — π0-corrected counts (helicity-resolved)

#include "pi0_corrected_counts.h"
#include "load_binning_scheme.h"  // Binning

#include <algorithm>
#include <cmath>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace fs = std::filesystem;

// ---------------- bin helpers (keep consistent with contamination stage) ----------------
using BinKey = std::tuple<int,int,int,int>; // (ix,iQ2,it,ip)

static constexpr int    N_PHI_BINS = 12;

static std::vector<std::pair<double,double>> uniqueRanges(const std::vector<Binning>& scheme, char which) {
    std::set<std::pair<double,double>> s;
    for (const auto& b : scheme) {
        if (which == 'x') s.emplace(b.xBmin, b.xBmax);
        else if (which == 'Q') s.emplace(b.Q2min, b.Q2max);
        else if (which == 't') s.emplace(b.tmin, b.tmax);
    }
    return std::vector<std::pair<double,double>>(s.begin(), s.end());
}

static inline std::string keyStr(const BinKey& k) {
    int a,b,c,d; std::tie(a,b,c,d)=k;
    std::ostringstream os; os<<"("<<a<<","<<b<<","<<c<<","<<d<<")";
    return os.str();
}

// ---------------- types for this stage ----------------
struct HelCounts { long long plus=0, minus=0; };          // raw DVCS counts
struct HelContam { double c_plus=0, c_plus_err=0, c_minus=0, c_minus_err=0; };
struct HelCorr   { double val_plus=0, err_plus=0, val_minus=0, err_minus=0; };

struct BinRecord {
    HelCounts N_data;    // from total_counts.json
    HelContam contam;    // from contamination_<period>.json
    HelCorr   N_corr;    // computed here
};

// ---------------- tiny JSON helpers ----------------
// Parse a signed/decimal/scientific number starting at or after ':' from 's' near position 'pos'.
static double parseNumberAfterColon(const std::string& s, size_t pos, double fallback = 0.0) {
    pos = s.find(':', pos);
    if (pos == std::string::npos) return fallback;
    size_t i = pos + 1;
    while (i < s.size() && std::isspace(static_cast<unsigned char>(s[i]))) ++i;
    size_t j = i;
    auto isnum=[&](char ch){
        return (std::isdigit(static_cast<unsigned char>(ch)) || ch=='-'||ch=='+'||ch=='.'||ch=='e'||ch=='E');
    };
    while (j < s.size() && isnum(s[j])) ++j;
    try { return std::stod(s.substr(i, j - i)); } catch (...) { return fallback; }
}

static long long parseIntAfterColon(const std::string& s, size_t pos, long long fallback = 0) {
    pos = s.find(':', pos);
    if (pos == std::string::npos) return fallback;
    size_t i = pos + 1;
    while (i < s.size() && std::isspace(static_cast<unsigned char>(s[i]))) ++i;
    size_t j = i;
    while (j < s.size() && (std::isdigit(static_cast<unsigned char>(s[j])) || s[j]=='-')) ++j;
    try { return std::stoll(s.substr(i, j - i)); } catch (...) { return fallback; }
}

// Return substring of the JSON object that corresponds to the value of the given key entry,
// when we already sit inside an object and have its body as 'obj'. Assumes "key": { ... } or "key": ....
static bool extractObjectForKey(const std::string& obj, const std::string& key, std::string& out) {
    size_t p = obj.find("\"" + key + "\"");
    if (p == std::string::npos) return false;
    p = obj.find(':', p); if (p == std::string::npos) return false;
    // If it's an object, capture balanced braces; otherwise capture up to next comma/brace
    size_t i = p + 1;
    while (i < obj.size() && std::isspace(static_cast<unsigned char>(obj[i]))) ++i;
    if (i < obj.size() && obj[i] == '{') {
        int depth = 0; size_t j = i;
        for (; j < obj.size(); ++j) {
            if (obj[j] == '{') depth++;
            else if (obj[j] == '}') { depth--; if (depth == 0) { ++j; break; } }
        }
        if (depth != 0) return false;
        out = obj.substr(i, j - i);
        return true;
    } else {
        // primitive value; take until comma or closing brace
        size_t j = i;
        while (j < obj.size() && obj[j] != ',' && obj[j] != '}') ++j;
        out = obj.substr(i, j - i);
        return true;
    }
}

// Iterate entries of a "bins": { "<key>": { ... }, ... } object and call fn(key,obj)
template <typename Fn>
static void iterateBins(const std::string& binsObj, Fn&& fn) {
    size_t kpos = 0;
    while (true) {
        size_t keyS = binsObj.find('"', kpos);
        if (keyS == std::string::npos) break;
        size_t keyE = binsObj.find('"', keyS+1);
        if (keyE == std::string::npos) break;
        std::string key = binsObj.substr(keyS+1, keyE-keyS-1); // e.g. "(ix,iQ2,it,ip)"

        size_t objS = binsObj.find('{', keyE);
        if (objS == std::string::npos) break;
        int depth=0; size_t j=objS;
        for (; j < binsObj.size(); ++j) {
            if (binsObj[j]=='{') depth++;
            else if (binsObj[j]=='}'){ depth--; if(!depth){ ++j; break; } }
        }
        if (depth != 0) break;
        std::string obj = binsObj.substr(objS, j-objS);

        fn(key, obj);
        kpos = j;
    }
}

// ---------------- readers ----------------
// total_counts.json reader (now per-period aware)
// Supports any of the following shapes:
//  A) { "bins": { "(ix,iQ2,it,ip)": { "N_data":{ "helicity":{"+1":..,"-1":..} } } } }
//  B) { "periods": { "<period>": { "bins": { ... } } } }
//  C) { "<period>": { "bins": { ... } }, "<other>": { ... } }   // rare, but be liberal
static std::map<BinKey, HelCounts> read_total_counts(const std::string& path,
                                                     const std::string& want_period)
{
    std::map<BinKey, HelCounts> out;

    std::ifstream ifs(path);
    if (!ifs) {
        std::cerr << "[pi0_corr][ERROR] Cannot open total counts JSON: " << path << std::endl;
        return out;
    }
    std::string s((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    ifs.close();

    auto findBalancedObject = [](const std::string& src, size_t fromBrace)->std::string {
        if (fromBrace == std::string::npos || fromBrace >= src.size() || src[fromBrace] != '{') return {};
        int depth = 0; size_t i = fromBrace;
        for (; i < src.size(); ++i) {
            if (src[i]=='{') depth++;
            else if (src[i]=='}'){ depth--; if (!depth){ ++i; break; } }
        }
        if (depth != 0) return {};
        return src.substr(fromBrace, i - fromBrace);
    };

    // Helper to parse a single "bins" object into out
    auto parseBinsObject = [&](const std::string& binsObj) {
        iterateBins(binsObj, [&](const std::string& key, const std::string& obj){
            int ix, iQ2, itb, ip;
            if (std::sscanf(key.c_str(), "(%d,%d,%d,%d)", &ix, &iQ2, &itb, &ip) != 4) return;

            long long Np = 0, Nm = 0;
            std::string NdataObj, helicityObj;

            if (extractObjectForKey(obj, "N_data", NdataObj)) {
                if (extractObjectForKey(NdataObj, "helicity", helicityObj)) {
                    size_t p1 = helicityObj.find("\"+1\"");
                    size_t p2 = helicityObj.find("\"-1\"");
                    if (p1 != std::string::npos) Np = parseIntAfterColon(helicityObj, p1, 0);
                    if (p2 != std::string::npos) Nm = parseIntAfterColon(helicityObj, p2, 0);
                }
            } else if (extractObjectForKey(obj, "helicity", helicityObj)) {
                size_t p1 = helicityObj.find("\"+1\"");
                size_t p2 = helicityObj.find("\"-1\"");
                if (p1 != std::string::npos) Np = parseIntAfterColon(helicityObj, p1, 0);
                if (p2 != std::string::npos) Nm = parseIntAfterColon(helicityObj, p2, 0);
            } else {
                size_t pP = obj.find("\"N_plus\"");
                size_t pM = obj.find("\"N_minus\"");
                if (pP != std::string::npos) Np = parseIntAfterColon(obj, pP, 0);
                if (pM != std::string::npos) Nm = parseIntAfterColon(obj, pM, 0);
            }

            out[BinKey(ix,iQ2,itb,ip)] = HelCounts{Np,Nm};
        });
    };

    // Case A: top-level "bins"
    {
        size_t pBins = s.find("\"bins\"");
        if (pBins != std::string::npos) {
            size_t brace = s.find('{', pBins);
            std::string binsObj = findBalancedObject(s, brace);
            if (!binsObj.empty()) {
                parseBinsObject(binsObj);
                if (!out.empty()) return out;
            }
        }
    }

    // Case B: "periods" -> "<period>" -> "bins"
    {
        size_t pPeriods = s.find("\"periods\"");
        if (pPeriods != std::string::npos) {
            size_t braceP = s.find('{', pPeriods);
            std::string periodsObj = findBalancedObject(s, braceP);
            if (!periodsObj.empty()) {
                // Find the requested period object
                std::string needle = "\"" + want_period + "\"";
                size_t pPer = periodsObj.find(needle);
                if (pPer != std::string::npos) {
                    size_t bracePer = periodsObj.find('{', pPer);
                    std::string perObj = findBalancedObject(periodsObj, bracePer);
                    if (!perObj.empty()) {
                        size_t pBins = perObj.find("\"bins\"");
                        if (pBins != std::string::npos) {
                            size_t brace = perObj.find('{', pBins);
                            std::string binsObj = findBalancedObject(perObj, brace);
                            if (!binsObj.empty()) {
                                parseBinsObject(binsObj);
                                if (!out.empty()) return out;
                            }
                        }
                    }
                }
            }
        }
    }

    // Case C: top-level "<period>" -> "bins"
    if (!want_period.empty()) {
        std::string needle = "\"" + want_period + "\"";
        size_t pTop = s.find(needle);
        if (pTop != std::string::npos) {
            size_t bracePer = s.find('{', pTop);
            std::string perObj = findBalancedObject(s, bracePer);
            if (!perObj.empty()) {
                size_t pBins = perObj.find("\"bins\"");
                if (pBins != std::string::npos) {
                    size_t brace = perObj.find('{', pBins);
                    std::string binsObj = findBalancedObject(perObj, brace);
                    if (!binsObj.empty()) {
                        parseBinsObject(binsObj);
                        if (!out.empty()) return out;
                    }
                }
            }
        }
    }

    std::cerr << "[pi0_corr][ERROR] Could not find a 'bins' object"
              << (want_period.empty() ? "" : (" for period '" + want_period + "'"))
              << " in " << path << "\n";
    return out;
}

// contamination_<period>.json: we rely on the format written by writeContaminationJson()
//   "contamination": {
//        "+1": { "value": <double>, "err": <double> },
//        "-1": { "value": <double>, "err": <double> }
//   }
static std::map<BinKey, HelContam> read_contamination(const std::string& path) {
    std::map<BinKey, HelContam> out;

    std::ifstream ifs(path);
    if (!ifs) {
        std::cerr << "[pi0_corr][ERROR] Cannot open contamination JSON: " << path << std::endl;
        return out;
    }
    std::string s((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    ifs.close();

    size_t pos = s.find("\"bins\"");
    if (pos == std::string::npos) { std::cerr << "[pi0_corr][ERROR] 'bins' not found in contamination JSON.\n"; return out; }
    pos = s.find('{', pos);
    if (pos == std::string::npos) return out;
    int depth = 0; size_t i = pos;
    for (; i < s.size(); ++i) { if (s[i]=='{') depth++; else if (s[i]=='}'){ depth--; if (!depth){ ++i; break; } } }
    std::string binsObj = s.substr(pos, i-pos);

    iterateBins(binsObj, [&](const std::string& key, const std::string& obj){
        int ix, iQ2, itb, ip;
        if (std::sscanf(key.c_str(), "(%d,%d,%d,%d)", &ix, &iQ2, &itb, &ip) != 4) return;

        std::string contObj;
        if (!extractObjectForKey(obj, "contamination", contObj)) return;

        double cP=0, eP=0, cM=0, eM=0;

        size_t pP = contObj.find("\"+1\"");
        if (pP != std::string::npos) {
            // read value/err inside "+1": { ... }
            size_t brace = contObj.find('{', pP);
            if (brace != std::string::npos) {
                int d=0; size_t j=brace; for(; j<contObj.size(); ++j){ if(contObj[j]=='{') d++; else if(contObj[j]=='}'){ d--; if(!d){ ++j; break; } } }
                std::string plusObj = contObj.substr(brace, j-brace);
                size_t pv = plusObj.find("\"value\"");
                size_t pe = plusObj.find("\"err\"");
                if (pv != std::string::npos) cP = parseNumberAfterColon(plusObj, pv, 0.0);
                if (pe != std::string::npos) eP = parseNumberAfterColon(plusObj, pe, 0.0);
            }
        }
        size_t pM = contObj.find("\"-1\"");
        if (pM != std::string::npos) {
            size_t brace = contObj.find('{', pM);
            if (brace != std::string::npos) {
                int d=0; size_t j=brace; for(; j<contObj.size(); ++j){ if(contObj[j]=='{') d++; else if(contObj[j]=='}'){ d--; if(!d){ ++j; break; } } }
                std::string minusObj = contObj.substr(brace, j-brace);
                size_t pv = minusObj.find("\"value\"");
                size_t pe = minusObj.find("\"err\"");
                if (pv != std::string::npos) cM = parseNumberAfterColon(minusObj, pv, 0.0);
                if (pe != std::string::npos) eM = parseNumberAfterColon(minusObj, pe, 0.0);
            }
        }

        out[BinKey(ix,iQ2,itb,ip)] = HelContam{cP,eP,cM,eM};
    });

    return out;
}

// ---------------- correction + writer ----------------
static HelCorr correct_one_helicity(long long N, double c, double c_err) {
    // N_corr = N * (1 - c)
    const double N_d   = static_cast<double>(std::max<long long>(0, N));
    const double one_mc = std::max(0.0, 1.0 - c);
    const double val = std::max(0.0, N_d * one_mc);

    // Var(N_corr) = (1 - c)^2 * Var(N) + (N)^2 * Var(c), with Var(N)=N (Poisson)
    const double var = (one_mc*one_mc) * N_d + (N_d*N_d) * (c_err*c_err);
    const double err = std::sqrt(std::max(0.0, var));

    HelCorr out;
    out.val_plus = val;  // NOTE: caller places into the right helicity field
    out.err_plus = err;
    return out;
}

static void write_corrected_json(
    const fs::path& out_path,
    const std::map<BinKey, BinRecord>& table,
    const std::vector<std::pair<double,double>>& xB_bins,
    const std::vector<std::pair<double,double>>& Q2_bins,
    const std::vector<std::pair<double,double>>& t_bins
) {
    std::ofstream ofs(out_path);
    if (!ofs) { std::cerr << "[pi0_corr][ERROR] Cannot open " << out_path << "\n"; return; }

    ofs << std::fixed << std::setprecision(8);
    ofs << "{\n";
    ofs << "  \"binning_meta\": {\"phi_bins\": " << N_PHI_BINS
        << ", \"xB_bins\": " << xB_bins.size()
        << ", \"Q2_bins\": " << Q2_bins.size()
        << ", \"t_bins\": "  << t_bins.size()  << "},\n";
    ofs << "  \"bins\": {\n";

    bool first=true;
    for (const auto& kv : table) {
        if (!first) ofs << ",\n";
        first=false;

        const auto& br = kv.second;
        const double tot_corr = br.N_corr.val_plus + br.N_corr.val_minus;
        const double tot_err  = std::sqrt(br.N_corr.err_plus*br.N_corr.err_plus
                                        + br.N_corr.err_minus*br.N_corr.err_minus);

        ofs << "    \"" << keyStr(kv.first) << "\": {";

        // original counts
        ofs << "\"N_data\":{\"helicity\":{\"+1\":" << br.N_data.plus
            << ",\"-1\":" << br.N_data.minus
            << "},\"total\":" << (br.N_data.plus + br.N_data.minus) << "}";

        // contamination inputs
        ofs << ",\"contamination\":{"
            << "\"+1\":{\"value\":" << br.contam.c_plus  << ",\"err\":" << br.contam.c_plus_err  << "},"
            << "\"-1\":{\"value\":" << br.contam.c_minus << ",\"err\":" << br.contam.c_minus_err << "}"
            << "}";

        // corrected counts
        ofs << ",\"N_corrected\":{"
            << "\"helicity\":{"
            << "\"+1\":{\"value\":" << br.N_corr.val_plus  << ",\"err\":" << br.N_corr.err_plus  << "},"
            << "\"-1\":{\"value\":" << br.N_corr.val_minus << ",\"err\":" << br.N_corr.err_minus << "}"
            << "},"
            << "\"total\":{\"value\":" << tot_corr << ",\"err\":" << tot_err << "}"
            << "}";

        ofs << "}";
    }

    ofs << "\n  }\n}\n";
    ofs.close();
    std::cout << "[pi0_corr] Wrote " << out_path.string() << std::endl;
}

static void write_combined_json(
    const fs::path& out_path,
    const std::map<std::string, std::map<BinKey, BinRecord>>& byPeriod,
    const std::vector<std::pair<double,double>>& xB_bins,
    const std::vector<std::pair<double,double>>& Q2_bins,
    const std::vector<std::pair<double,double>>& t_bins
) {
    std::ofstream ofs(out_path);
    if (!ofs) { std::cerr << "[pi0_corr][ERROR] Cannot open combined output " << out_path << "\n"; return; }

    ofs << std::fixed << std::setprecision(8);
    ofs << "{\n";
    ofs << "  \"binning_meta\": {\"phi_bins\": " << N_PHI_BINS
        << ", \"xB_bins\": " << xB_bins.size()
        << ", \"Q2_bins\": " << Q2_bins.size()
        << ", \"t_bins\": "  << t_bins.size()  << "},\n";
    ofs << "  \"periods\": {\n";

    bool firstP = true;
    for (const auto& pkv : byPeriod) {
        if (!firstP) ofs << ",\n";
        firstP = false;

        ofs << "    \"" << pkv.first << "\": {\n";
        ofs << "      \"bins\": {\n";

        bool firstB = true;
        for (const auto& kv : pkv.second) {
            if (!firstB) ofs << ",\n";
            firstB = false;

            const auto& br = kv.second;
            const double tot_corr = br.N_corr.val_plus + br.N_corr.val_minus;
            const double tot_err  = std::sqrt(br.N_corr.err_plus*br.N_corr.err_plus
                                            + br.N_corr.err_minus*br.N_corr.err_minus);

            ofs << "        \"" << keyStr(kv.first) << "\": {";
            ofs << "\"N_data\":{\"helicity\":{\"+1\":" << br.N_data.plus
                << ",\"-1\":" << br.N_data.minus << "},\"total\":" << (br.N_data.plus + br.N_data.minus) << "}";
            ofs << ",\"contamination\":{"
                << "\"+1\":{\"value\":" << br.contam.c_plus  << ",\"err\":" << br.contam.c_plus_err  << "},"
                << "\"-1\":{\"value\":" << br.contam.c_minus << ",\"err\":" << br.contam.c_minus_err << "}"
                << "}";
            ofs << ",\"N_corrected\":{"
                << "\"helicity\":{"
                << "\"+1\":{\"value\":" << br.N_corr.val_plus  << ",\"err\":" << br.N_corr.err_plus  << "},"
                << "\"-1\":{\"value\":" << br.N_corr.val_minus << ",\"err\":" << br.N_corr.err_minus << "}"
                << "},"
                << "\"total\":{\"value\":" << tot_corr << ",\"err\":" << tot_err << "}"
                << "}";
            ofs << "}";
        }

        ofs << "\n      }\n"; // bins
        ofs << "    }";
    }

    ofs << "\n  }\n}\n";
    ofs.close();
    std::cout << "[pi0_corr] Wrote combined " << out_path.string() << std::endl;
}

// ---------------- main entry ----------------
void compute_pi0_corrected_counts(
    const std::vector<std::string>& periods,
    const std::vector<Binning>& binning_scheme,
    const std::string& total_counts_json,
    const std::string& contamination_dir,
    const std::string& out_root_dir
) {
    // binning meta (for JSON headers only)
    const auto xB_bins = uniqueRanges(binning_scheme, 'x');
    const auto Q2_bins = uniqueRanges(binning_scheme, 'Q');
    const auto t_bins  = uniqueRanges(binning_scheme, 't');
    if (xB_bins.empty() || Q2_bins.empty() || t_bins.empty()) {
        std::cerr << "[pi0_corr][ERROR] Missing binning ranges." << std::endl;
        return;
    }

    if (!fs::exists(total_counts_json)) {
        std::cerr << "[pi0_corr][ERROR] total_counts.json not found: " << total_counts_json << std::endl;
        return;
    }
    if (!fs::exists(contamination_dir)) {
        std::cerr << "[pi0_corr][ERROR] contamination dir not found: " << contamination_dir << std::endl;
        return;
    }

    // Prepare output dirs
    const fs::path out_root(out_root_dir);
    const fs::path out_json_dir = out_root / "jsons";
    std::error_code ec;
    fs::create_directories(out_json_dir, ec);

    // Keep all per-period tables to write a combined file
    std::map<std::string, std::map<BinKey, BinRecord>> allPeriods;

    // For each period: read contamination JSON and produce corrected counts
    for (const auto& period : periods) {
        const fs::path contam_path = fs::path(contamination_dir) / ("contamination_" + period + ".json");
        if (!fs::exists(contam_path)) {
            std::cerr << "[pi0_corr][WARN] Missing contamination JSON for " << period
                      << " (" << contam_path.string() << "). Skipping period." << std::endl;
            continue;
        }
        auto contam_map = read_contamination(contam_path.string());
        if (contam_map.empty()) {
            std::cerr << "[pi0_corr][WARN] No bins in contamination JSON for " << period << ". Skipping period." << std::endl;
            continue;
        }

        // NEW: read total counts *for this period*
        const auto total_counts_map = read_total_counts(total_counts_json, period);
        if (total_counts_map.empty()) {
            std::cerr << "[pi0_corr][WARN] No DVCS counts found in total_counts.json for period '"
                      << period << "'. Skipping period." << std::endl;
            continue;
        }

        std::map<BinKey, BinRecord> table;

        for (const auto& kv : total_counts_map) {
            const BinKey& key = kv.first;
            BinRecord rec;
            rec.N_data = kv.second;

            auto itC = contam_map.find(key);
            if (itC != contam_map.end())
                rec.contam = itC->second;

            // +1 helicity
            {
                const auto tmp = correct_one_helicity(rec.N_data.plus, rec.contam.c_plus, rec.contam.c_plus_err);
                rec.N_corr.val_plus = tmp.val_plus;
                rec.N_corr.err_plus = tmp.err_plus;
            }
            // -1 helicity
            {
                const auto tmp = correct_one_helicity(rec.N_data.minus, rec.contam.c_minus, rec.contam.c_minus_err);
                rec.N_corr.val_minus = tmp.val_plus;
                rec.N_corr.err_minus = tmp.err_plus;
            }

            table[key] = rec;
        }

        // Write & stash
        const fs::path out_period = out_json_dir / ("pi0_corrected_counts_" + period + ".json");
        write_corrected_json(out_period, table, xB_bins, Q2_bins, t_bins);
        allPeriods[period] = std::move(table);
    }

    // Combined file
    const fs::path out_combined = out_json_dir / "pi0_corrected_counts_all_periods.json";
    write_combined_json(out_combined, allPeriods, xB_bins, Q2_bins, t_bins);
}

and finally a bsa.cpp:

// bsa.cpp
#include "bsa.h"

#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TAxis.h>
#include <TStyle.h>
#include <TH1.h>
#include <TF1.h>
#include <TLatex.h>
#include <TFitResult.h>
#include <TPad.h>
#include <TGaxis.h>
#include <TTree.h>

#include <Math/Factory.h>
#include <Math/Minimizer.h>
#include <Math/Functor.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace {

// =================== global style ===================
struct StyleInit {
    StyleInit() {
        gStyle->SetOptTitle(0);
        gStyle->SetOptStat(0);
        gStyle->SetFrameLineWidth(2);
        gStyle->SetLineWidth(2);
        gStyle->SetPadTickX(1);
        gStyle->SetPadTickY(1);
        gStyle->SetLegendBorderSize(1);

        const int rf = 42; // Helvetica
        gStyle->SetTitleFont(rf, "XYZ");
        gStyle->SetLabelFont(rf, "XYZ");
        gStyle->SetTextFont(rf);
    }
} _style_bootstrap;

constexpr int    N_PHI_BINS = 12;
constexpr double TWO_PI     = 2.0 * M_PI;

// ---- BSA fit stabilization knobs (adapted to single-cos denominator) ----
constexpr double EPS_DEN_FLOOR = 1e-2;   // demand D(φ) ≥ this over [0,2π)
constexpr double LAMBDA_DEN    = 1e6;    // penalty multiplier for denominator floor
constexpr double EPS_DEN_EVAL  = 1e-6;   // evaluation clamp for denominator

constexpr double A_MAX_AMP     = 0.999;  // soft "box" for |B| (and optional |A|,|C|)
constexpr double LAMBDA_AMP    = 1e5;    // penalty multiplier for amplitude box

using BinKey = std::tuple<int,int,int,int>; // (ix,iQ2,it,ip)
struct HelCounts { long long plus=0, minus=0; };

static inline std::string toLower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c){ return std::tolower(c); });
    return s;
}
static std::string periodToRunTagKey(const std::string& period) {
    auto pos = period.find('_');
    if (pos == std::string::npos || pos + 1 >= period.size()) return toLower(period);
    return toLower(period.substr(pos + 1));
}
static std::vector<std::pair<double,double>> uniqueRanges(const std::vector<Binning>& scheme, char which) {
    std::set<std::pair<double,double>> s;
    for (const auto& b : scheme) {
        if (which == 'x') s.emplace(b.xBmin, b.xBmax);
        else if (which == 'Q') s.emplace(b.Q2min, b.Q2max);
        else if (which == 't') s.emplace(b.tmin, b.tmax);
    }
    return std::vector<std::pair<double,double>>(s.begin(), s.end());
}
static inline std::vector<double> phiCentersRad() {
    std::vector<double> v(N_PHI_BINS);
    const double step = TWO_PI / double(N_PHI_BINS);
    for (int i=0;i<N_PHI_BINS;++i) v[i] = (i+0.5)*step;
    return v;
}
static inline std::vector<double> phiCentersDeg() {
    std::vector<double> d(N_PHI_BINS);
    const double step = 360.0/double(N_PHI_BINS);
    for (int i=0;i<N_PHI_BINS;++i) d[i] = (i+0.5)*step;
    return d;
}
static inline int findIndex(const std::pair<double,double>& range,
                            const std::vector<std::pair<double,double>>& ranges) {
    for (int i=0;i<(int)ranges.size();++i) if (ranges[i]==range) return i;
    return -1;
}

// For a given xB range, list only the Q² and |t| ranges present in the CSV/binning.
static void uniqueQT_for_xB(
    const std::vector<Binning>& scheme,
    const std::pair<double,double>& xBrange,
    std::vector<std::pair<double,double>>& Q2_list,
    std::vector<std::pair<double,double>>& t_list
) {
    std::set<std::pair<double,double>> qs, ts;
    for (const auto& b : scheme) {
        if (std::make_pair(b.xBmin,b.xBmax) == xBrange) {
            qs.emplace(b.Q2min,b.Q2max);
            ts.emplace(b.tmin,b.tmax);
        }
    }
    Q2_list.assign(qs.begin(), qs.end());
    t_list.assign(ts.begin(), ts.end());
}

// ------------ I/O helpers: total_counts.json ------------
using GroupCounts = std::map<std::string, std::map<BinKey, HelCounts>>;

static bool parse_tuple_key(const std::string& s, BinKey& out) {
    int ix,iQ,it,ip;
    if (std::sscanf(s.c_str(),"(%d,%d,%d,%d)",&ix,&iQ,&it,&ip)!=4) return false;
    out = BinKey(ix,iQ,it,ip);
    return true;
}

// Minimal parser (shape written by your total_counts code)
static bool load_total_counts(const std::string& path, GroupCounts& outGroups) {
    std::ifstream ifs(path);
    if (!ifs) { std::cerr<<"[bsa][ERROR] Cannot open total_counts JSON: "<<path<<"\n"; return false; }
    std::string s((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());

    size_t gpos = s.find("\"groups\"");
    if (gpos==std::string::npos) { std::cerr<<"[bsa][ERROR] 'groups' not found in total_counts.\n"; return false; }
    size_t brace = s.find('{', gpos); if (brace==std::string::npos) return false;
    int d=0; size_t i=brace; for (; i<s.size(); ++i){ if(s[i]=='{') d++; else if(s[i]=='}'){ d--; if(!d){ ++i; break; } } }
    std::string groupsObj = s.substr(brace, i-brace);

    size_t kpos=0;
    while (true) {
        size_t q1 = groupsObj.find('"', kpos); if (q1==std::string::npos) break;
        size_t q2 = groupsObj.find('"', q1+1); if (q2==std::string::npos) break;
        std::string gname = groupsObj.substr(q1+1, q2-q1-1);

        size_t objS = groupsObj.find('{', q2); if (objS==std::string::npos) break;
        int d2=0; size_t j=objS;
        for (; j<groupsObj.size(); ++j){ if(groupsObj[j]=='{') d2++; else if(groupsObj[j]=='}'){ d2--; if(!d2){ ++j; break; } } }
        std::string binsObj = groupsObj.substr(objS, j-objS);

        std::map<BinKey, HelCounts> gmap;
        size_t bpos=0;
        while (true) {
            size_t bk1 = binsObj.find('"', bpos); if (bk1==std::string::npos) break;
            size_t bk2 = binsObj.find('"', bk1+1); if (bk2==std::string::npos) break;
            std::string key = binsObj.substr(bk1+1, bk2-bk1-1);
            BinKey bk;
            if (!parse_tuple_key(key, bk)) { bpos=bk2+1; continue; }

            size_t valS = binsObj.find('{', bk2); if (valS==std::string::npos) break;
            int d3=0; size_t jj=valS;
            for (; jj<binsObj.size(); ++jj){ if(binsObj[jj]=='{') d3++; else if(binsObj[jj]=='}'){ d3--; if(!d3){ ++jj; break; } } }
            std::string obj = binsObj.substr(valS, jj-valS);

            auto findLL = [&](const char* pat)->long long {
                size_t p = obj.find(pat); if (p==std::string::npos) return 0;
                p = obj.find(':', p); if (p==std::string::npos) return 0;
                size_t a=p+1; while (a<obj.size() && isspace((unsigned char)obj[a])) ++a;
                size_t b=a; while (b<obj.size() && (isdigit((unsigned char)obj[b])||obj[b]=='-')) ++b;
                try { return std::stoll(obj.substr(a,b-a)); } catch(...) { return 0; }
            };
            HelCounts hc;
            hc.plus  = findLL("\"+1\"");
            hc.minus = findLL("\"-1\"");

            gmap[bk]=hc;
            bpos=jj;
        }

        outGroups[gname]=std::move(gmap);
        kpos=j;
    }
    return !outGroups.empty();
}

// ------------ I/O helpers: contamination_<period>.json ------------
struct ContamBin { double c_plus=0, c_plus_err=0, c_minus=0, c_minus_err=0; };
using ContamMap = std::map<BinKey, ContamBin>;

static bool load_contam_for_period(const std::string& path, ContamMap& out) {
    std::ifstream ifs(path);
    if (!ifs) return false;
    std::string s((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    size_t pos = s.find("\"bins\""); if (pos==std::string::npos) return false;
    size_t br = s.find('{', pos); if (br==std::string::npos) return false;
    int d=0; size_t i=br; for (; i<s.size(); ++i){ if(s[i]=='{') d++; else if(s[i]=='}'){ d--; if(!d){ ++i; break; } } }
    std::string binsObj = s.substr(br, i-br);

    size_t kpos=0;
    while (true) {
        size_t q1 = binsObj.find('"', kpos); if (q1==std::string::npos) break;
        size_t q2 = binsObj.find('"', q1+1); if (q2==std::string::npos) break;
        std::string key = binsObj.substr(q1+1, q2-q1-1);
        BinKey bk; if (!parse_tuple_key(key, bk)) { kpos=q2+1; continue; }

        size_t objS = binsObj.find('{', q2); if (objS==std::string::npos) break;
        int d2=0; size_t j=objS; for (; j<binsObj.size(); ++j){ if(binsObj[j]=='{') d2++; else if(binsObj[j]=='}'){ d2--; if(!d2){ ++j; break; } } }
        std::string obj = binsObj.substr(objS, j-objS);

        auto findD = [&](const char* pat)->double{
            size_t p=obj.find(pat); if (p==std::string::npos) return 0.0;
            p=obj.find(':',p); if (p==std::string::npos) return 0.0;
            size_t a=p+1; while (a<obj.size() && isspace((unsigned char)obj[a])) ++a;
            size_t b=a; while (b<obj.size() && (isdigit((unsigned char)obj[b])||obj[b]=='-'||obj[b]=='.'||obj[b]=='e'||obj[b]=='E'||obj[b]=='+')) ++b;
            try { return std::stod(obj.substr(a,b-a)); } catch(...) { return 0.0; }
        };
        ContamBin cb;
        cb.c_plus      = findD("\"contamination\":{\"+1\":{\"value\"");
        cb.c_plus_err  = findD("\"contamination\":{\"+1\":{\"err\"");
        cb.c_minus     = findD("\"contamination\":{\"-1\":{\"value\"");
        cb.c_minus_err = findD("\"contamination\":{\"-1\":{\"err\"");
        out[bk]=cb;
        kpos=j;
    }
    return true;
}

// ------------ polarization (from DVCS tree) ------------
struct PolStats { double P=1.0; double P_se=0.0; int n=0; };

struct BranchPol {
    int helicity=0;
    double t1=0, open_angle_ep2=0, pTmiss=0;
    double x=0, Q2=0, phi2=0, Delta_phi=0;
    double beam_pol=1.0;
    int detector1=0, detector2=0;
    bool hasH=false, hasCuts=false, hasBins=false, hasPhi2=false, hasDp=false, hasPol=false, hasTopo=false;
    void bind(TTree* t){
        if (!t) return;
        auto bindI=[&](const char*n,int*a,bool&f){ if(t->GetBranch(n)){ t->SetBranchAddress(n,a); f=true; } };
        auto bindD=[&](const char*n,double*a,bool&f){ if(t->GetBranch(n)){ t->SetBranchAddress(n,a); f=true; } };
        bindI("helicity",&helicity,hasH);
        bindD("t1",&t1,hasCuts); bindD("open_angle_ep2",&open_angle_ep2,hasCuts); bindD("pTmiss",&pTmiss,hasCuts);
        bindD("x",&x,hasBins); bindD("Q2",&Q2,hasBins); bindD("phi2",&phi2,hasPhi2); bindD("Delta_phi",&Delta_phi,hasDp);
        bindD("beam_pol",&beam_pol,hasPol);
        bindI("detector1",&detector1,hasTopo); bindI("detector2",&detector2,hasTopo);
    }
    double phi() const { return hasPhi2 ? phi2 : (hasDp ? Delta_phi : std::numeric_limits<double>::quiet_NaN()); }
};

static inline bool passesTopology_simple(int d1, int d2, const std::vector<std::string>& tops){
    for (const auto& t : tops){
        if (t=="(FD,FD)" && d1==1 && d2==1) return true;
        if (t=="(CD,FD)" && d1==2 && d2==1) return true;
        if (t=="(CD,FT)" && d1==2 && d2==0) return true;
    }
    return false;
}
static inline bool applyKinematicCuts_simple(double t1,double oa,double pT){ return !(oa<=5.0 || (-t1)>1.0 || pT>0.20); }

static std::vector<PolStats> compute_bin_polarization(
    TTree* t, const std::vector<std::pair<double,double>>& xB_bins,
    const std::vector<std::pair<double,double>>& Q2_bins,
    const std::vector<std::pair<double,double>>& t_bins,
    const std::vector<std::string>& tops)
{
    std::vector<PolStats> Pbin(xB_bins.size()*Q2_bins.size()*t_bins.size()*N_PHI_BINS);
    auto idx=[&](int ix,int iQ,int it,int ip){ return (((ix* (int)Q2_bins.size() + iQ)*(int)t_bins.size()+it)*N_PHI_BINS + ip); };

    BranchPol b; b.bind(t);
    if (!(b.hasH && b.hasCuts && b.hasBins && (b.hasPhi2||b.hasDp) && b.hasPol && b.hasTopo)) return Pbin;

    const Long64_t nent = t->GetEntries();
    for (Long64_t i=0;i<nent;++i){
        t->GetEntry(i);
        if (!applyKinematicCuts_simple(b.t1,b.open_angle_ep2,b.pTmiss)) continue;
        if (b.helicity!=+1 && b.helicity!=-1) continue;
        if (!passesTopology_simple(b.detector1,b.detector2,tops)) continue;

        double xB=b.x, Q2=b.Q2, tt=fabs(b.t1), phi=b.phi();
        if (!std::isfinite(xB)||!std::isfinite(Q2)||!std::isfinite(tt)||!std::isfinite(phi)) continue;

        auto findBin=[&](double v,const std::vector<std::pair<double,double>>& ranges)->int{
            for (int k=0;k<(int)ranges.size();++k) if (v>=ranges[k].first && v<ranges[k].second) return k;
            return -1;
        };
        int ix=findBin(xB,xB_bins), iQ=findBin(Q2,Q2_bins), it=findBin(tt,t_bins);
        if (ix<0||iQ<0||it<0) continue;

        double width = TWO_PI/double(N_PHI_BINS);
        double w = std::fmod(phi,TWO_PI); if (w<0) w+=TWO_PI; if (w>=TWO_PI) w = std::nextafter(TWO_PI,0.0);
        int ip = std::min(std::max(int(std::floor(w/width)),0), N_PHI_BINS-1);

        PolStats& ps = Pbin[idx(ix,iQ,it,ip)];
        ps.P += (b.beam_pol - ps.P)/(ps.n+1); // online mean
        ps.n++;
    }

    for (auto& ps : Pbin) ps.P_se = 0.0;
    return Pbin;
}

// ------------ math helpers for stabilizer penalties ------------
static inline double Dmin_single(double B) {
    // For D(φ)=1 + B cosφ, min over φ is 1 - |B|
    return 1.0 - std::fabs(B);
}
static inline double overBox(double A){ double v = std::fabs(A) - A_MAX_AMP; return (v>0.0)? v*v : 0.0; }

// ------------ BSA computation core ------------
struct BSApt { double phi=0.0; double bsa=0.0; double err=0.0; bool valid=false; };
struct FitRes { double A=0, Aerr=0, B1=0, B1err=0, B2=0, B2err=0, C=0, Cerr=0, chi2=0; int ndf=0; int status=0; };
struct CellResult {
    std::vector<BSApt> points; FitRes fit;
    double P_used=1.0; bool P_per_bin=false; double P_period_avg=1.0;
};

// Custom χ² + penalties minimization (Minuit2) for A, B(=B1), C
static FitRes fit_cell(const std::vector<BSApt>& pts){
    std::vector<double> x, y, ey;
    x.reserve(pts.size()); y.reserve(pts.size()); ey.reserve(pts.size());
    for (const auto& p : pts) {
        if (!p.valid) continue;
        x.push_back(p.phi);                 // radians
        y.push_back(p.bsa);
        ey.push_back(std::max(p.err, 1e-6));
    }
    FitRes fr; 
    const int n = (int)x.size();
    if (n < 4) { fr.status = 1; fr.ndf = 0; return fr; }

    // χ² + stabilizers
    auto chi2pen = [&](const double *par){
        const double A  = par[0];
        const double B1 = par[1]; // single cosine coefficient
        const double C  = par[2];

        double chi2 = 0.0;
        for (int i=0;i<n;++i){
            const double phi = x[i];
            double denom = 1.0 + B1*std::cos(phi);
            if (denom < EPS_DEN_EVAL) denom = EPS_DEN_EVAL; // evaluation clamp
            const double yhat = C + (A*std::sin(phi))/denom;
            const double pull = (y[i] - yhat)/ey[i];
            chi2 += pull*pull;
        }

        // (1) Denominator soft barrier over full φ-domain: min = 1 - |B|
        const double Dmin = Dmin_single(B1);
        double pen_den = 0.0;
        if (Dmin < EPS_DEN_FLOOR) {
            const double deficit = EPS_DEN_FLOOR - Dmin;
            pen_den = LAMBDA_DEN * deficit * deficit;
        }

        // (2) Amplitude box on |B|
        double pen_amp = overBox(B1);
        // Optionally also constrain |A|,|C| softly:
        // pen_amp += 0.25*overBox(A); pen_amp += 0.25*overBox(C);

        return chi2 + pen_den + LAMBDA_AMP*pen_amp;
    };

    // Minimizer
    std::unique_ptr<ROOT::Math::Minimizer> min(
        ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad"));
    min->SetMaxFunctionCalls(10000);
    min->SetMaxIterations(10000);
    min->SetTolerance(1e-6);
    min->SetPrintLevel(0);

    ROOT::Math::Functor fcn(chi2pen, 3);
    min->SetFunction(fcn);

    // Parameters: A, B (as B1), C
    min->SetLimitedVariable(0, "A",  0.10, 0.02, -1.0,  1.0);
    min->SetLimitedVariable(1, "B1", 0.00, 0.02, -1.0,  1.0);
    min->SetLimitedVariable(2, "C",  0.00, 0.02, -0.5,  0.5);

    const bool ok = min->Minimize();
    fr.status = ok ? 0 : 1;

    const double *par = min->X();
    const double *err = min->Errors();

    fr.A  = par[0]; fr.Aerr  = err[0];
    fr.B1 = par[1]; fr.B1err = err[1];
    fr.B2 = 0.0;    fr.B2err = 0.0; // single-cos fit => no cos2φ term
    fr.C  = par[2]; fr.Cerr  = err[2];

    fr.chi2 = chi2pen(par);
    fr.ndf  = std::max(0, n - 3);
    return fr;
}

// JSON writers (schema kept: B1 used, B2 always 0)
static void write_period_bsa_json(
    const std::string& out_path,
    int nPhi,
    const std::vector<std::pair<double,double>>& xB_bins,
    const std::vector<std::pair<double,double>>& Q2_bins,
    const std::vector<std::pair<double,double>>& t_bins,
    const std::map<std::tuple<int,int,int>, CellResult>& cells)
{
    std::ofstream ofs(out_path);
    if (!ofs) { std::cerr<<"[bsa][ERROR] Cannot open "<<out_path<<"\n"; return; }
    ofs<<std::fixed<<std::setprecision(8);
    ofs<<"{\n";
    ofs<<"  \"binning_meta\": {\"phi_bins\": "<<nPhi<<", \"xB_bins\": "<<xB_bins.size()<<", \"Q2_bins\": "<<Q2_bins.size()<<", \"t_bins\": "<<t_bins.size()<<"},\n";
    ofs<<"  \"bins\": {\n";
    bool first=true;
    for (const auto& kv : cells){
        if (!first) ofs<<",\n"; first=false;
        int ix,iQ,it; std::tie(ix,iQ,it)=kv.first;
        const auto& cr = kv.second;
        ofs<<"    \"("<<ix<<","<<iQ<<","<<it<<")\": {\n";
        ofs<<"      \"phi\": [";
        for (size_t i=0;i<cr.points.size();++i){ if(i) ofs<<","; ofs<<cr.points[i].phi; } ofs<<"],\n";
        ofs<<"      \"bsa\": [";
        for (size_t i=0;i<cr.points.size();++i){ if(i) ofs<<","; ofs<<cr.points[i].bsa; } ofs<<"],\n";
        ofs<<"      \"bsa_err\": [";
        for (size_t i=0;i<cr.points.size();++i){ if(i) ofs<<","; ofs<<cr.points[i].err; } ofs<<"],\n";
        ofs<<"      \"fit\": {"
              "\"A\": {\"val\": "<<cr.fit.A<<", \"err\": "<<cr.fit.Aerr<<"}, "
              "\"B1\": {\"val\": "<<cr.fit.B1<<", \"err\": "<<cr.fit.B1err<<"}, "
              "\"B2\": {\"val\": 0.0, \"err\": 0.0}, "
              "\"C\": {\"val\": "<<cr.fit.C<<", \"err\": "<<cr.fit.Cerr<<"}, "
              "\"chi2\": "<<cr.fit.chi2<<", \"ndf\": "<<cr.fit.ndf<<", \"status\": "<<cr.fit.status<<"},\n";
        ofs<<"      \"polarization\": {\"per_bin\": "<<(cr.P_per_bin?"true":"false")<<", \"P_used\": "<<cr.P_used<<"}\n";
        ofs<<"    }";
    }
    ofs<<"\n  }\n}\n";
}

static void write_all_periods_json(
    const std::string& out_path,
    const std::map<std::string, std::map<std::tuple<int,int,int>, CellResult>>& perPeriod,
    int nPhi,
    const std::vector<std::pair<double,double>>& xB_bins,
    const std::vector<std::pair<double,double>>& Q2_bins,
    const std::vector<std::pair<double,double>>& t_bins)
{
    std::ofstream ofs(out_path);
    if (!ofs) { std::cerr<<"[bsa][ERROR] Cannot open "<<out_path<<"\n"; return; }
    ofs<<std::fixed<<std::setprecision(8);
    ofs<<"{\n";
    ofs<<"  \"binning_meta\": {\"phi_bins\": "<<nPhi<<", \"xB_bins\": "<<xB_bins.size()<<", \"Q2_bins\": "<<Q2_bins.size()<<", \"t_bins\": "<<t_bins.size()<<"},\n";
    ofs<<"  \"periods\": {\n";
    bool firstP=true;
    for (const auto& pkv : perPeriod){
        if (!firstP) ofs<<",\n"; firstP=false;
        ofs<<"    \""<<pkv.first<<"\": {\n      \"bins\": {\n";
        bool firstB=true;
        for (const auto& kv : pkv.second){
            if (!firstB) ofs<<",\n"; firstB=false;
            int ix,iQ,it; std::tie(ix,iQ,it)=kv.first;
            const auto& cr = kv.second;
            ofs<<"        \"("<<ix<<","<<iQ<<","<<it<<")\": {";
            ofs<<"\"phi\":["; for (size_t i=0;i<cr.points.size();++i){ if(i)ofs<<","; ofs<<cr.points[i].phi; } ofs<<"],";
            ofs<<"\"bsa\":["; for (size_t i=0;i<cr.points.size();++i){ if(i)ofs<<","; ofs<<cr.points[i].bsa; } ofs<<"],";
            ofs<<"\"bsa_err\":["; for (size_t i=0;i<cr.points.size();++i){ if(i)ofs<<","; ofs<<cr.points[i].err; } ofs<<"],";
            ofs<<"\"fit\":{"
                  "\"A\":{\"val\":"<<cr.fit.A<<",\"err\":"<<cr.fit.Aerr<<"},"
                  "\"B1\":{\"val\":"<<cr.fit.B1<<",\"err\":"<<cr.fit.B1err<<"},"
                  "\"B2\":{\"val\":0.0,\"err\":0.0},"
                  "\"C\":{\"val\":"<<cr.fit.C<<",\"err\":"<<cr.fit.Cerr<<"},"
                  "\"chi2\":"<<cr.fit.chi2<<",\"ndf\":"<<cr.fit.ndf<<",\"status\":"<<cr.fit.status<<"},";
            ofs<<"\"polarization\":{\"per_bin\":"<<(cr.P_per_bin?"true":"false")<<",\"P_used\":"<<cr.P_used<<"}";
            ofs<<"}";
        }
        ofs<<"\n      }\n    }";
    }
    ofs<<"\n  }\n}\n";
}

// ------------ plotting (dynamic grid + clean title + red fit curve) ------------
static void drawDegreeTicks(double xmin, double ymin, double xmax, double labelSize){
    // Draw 0,90,180,270,360 degree ticks with labels on top of the frame axis
    TGaxis* ax = new TGaxis(xmin, ymin, xmax, ymin, 0.0, 360.0, 4, "");
    ax->SetLabelFont(42);
    ax->SetLabelSize(labelSize);
    ax->SetLabelOffset(0.012);
    ax->SetTitle("");
    ax->SetTickSize(0.02);
    ax->Draw();
}

static void plot_cells_for_period(
    const std::string& period,
    const std::vector<Binning>& binning_scheme,
    const std::vector<std::pair<double,double>>& xB_bins,
    const std::vector<std::pair<double,double>>& Q2_bins,
    const std::vector<std::pair<double,double>>& t_bins,
    const std::map<std::tuple<int,int,int>, CellResult>& cells,
    const std::string& out_dir_plots)
{
    using std::filesystem::create_directories;
    std::error_code ec;
    create_directories(out_dir_plots, ec);

    static const auto PHI_DEG = phiCentersDeg();

    for (int ix = 0; ix < (int)xB_bins.size(); ++ix) {
        const auto xb = xB_bins[ix];

        // Only the Q² and |t| bins that exist in this xB slice
        std::vector<std::pair<double,double>> Q2_slice, t_slice;
        uniqueQT_for_xB(binning_scheme, xb, Q2_slice, t_slice);
        if (Q2_slice.empty() || t_slice.empty()) continue;

        const int nrows = (int)t_slice.size();
        const int ncols = (int)Q2_slice.size();

        // Separate title pad (smaller now)
        const int W = 280*ncols + 160;
        const int H = 240*nrows + 170;

        std::ostringstream cname; cname<<"c_bsa_"<<period<<"_xB"<<ix;
        TCanvas* c = new TCanvas(cname.str().c_str(), cname.str().c_str(), W, H);

        TPad* pTop  = new TPad("pTop","pTop", 0.0, 0.915, 1.0, 1.0);
        pTop->SetFillStyle(0); pTop->SetBorderSize(0); pTop->Draw();

        TPad* pGrid = new TPad("pGrid","pGrid", 0.0, 0.00, 1.0, 0.915);
        pGrid->SetFillStyle(0); pGrid->SetBorderSize(0); pGrid->Draw();
        pGrid->cd();
        pGrid->Divide(ncols, nrows, 0.0001, 0.0001);

        // Title (compact and safe)
        pTop->cd();
        TLatex head;
        head.SetNDC(); head.SetTextAlign(22);
        head.SetTextFont(42);
        head.SetTextSize(0.36); // smaller than before
        std::ostringstream tit;
        tit << Form("Beam-Spin Asymmetry  %s   x_{B} #in [%.2g, %.2g]",
                    period.c_str(), xb.first, xb.second);
        head.DrawLatex(0.5, 0.55, tit.str().c_str());

        // Panels
        for (int r = 0; r < nrows; ++r) {
            const int it_global = findIndex(t_slice[r], t_bins);
            if (it_global < 0) continue;

            for (int ccol = 0; ccol < ncols; ++ccol) {
                const int iQ_global = findIndex(Q2_slice[ccol], Q2_bins);
                if (iQ_global < 0) continue;

                pGrid->cd(r*ncols + ccol + 1);
                gPad->SetGrid(1,1);
                gPad->SetTopMargin(0.08);
                gPad->SetBottomMargin(0.18);
                gPad->SetLeftMargin(0.15);
                gPad->SetRightMargin(0.10);

                // frame
                TH1* frame = gPad->DrawFrame(0.0, -1.05, 360.0, 1.05);
                TAxis* ax = frame->GetXaxis();
                TAxis* ay = frame->GetYaxis();

                // We will overlay a TGaxis for 0,90,180,270,360, so hide the default labels
                ax->SetLabelSize(0.0001);
                ax->SetTitle("#phi (deg)");
                ay->SetTitle("A_{LU}");
                ax->CenterTitle(); ay->CenterTitle();
                ax->SetNdivisions(505);

                ax->SetTitleSize(0.060);
                ay->SetTitleSize(0.060);
                ay->SetLabelSize(0.048);

                ax->SetTitleOffset(1.25);
                ay->SetTitleOffset(1.35);

                // overlay the custom degree ticks
                drawDegreeTicks(gPad->GetUxmin(), gPad->GetUymin(), gPad->GetUxmax(), 0.048);

                auto itCell = cells.find(std::make_tuple(ix, iQ_global, it_global));
                if (itCell == cells.end()) continue;
                const auto& cr = itCell->second;

                // Data points
                std::vector<double> x, y, ey;
                x.reserve(N_PHI_BINS); y.reserve(N_PHI_BINS); ey.reserve(N_PHI_BINS);
                for (int ip=0; ip<N_PHI_BINS; ++ip) {
                    const auto& p = cr.points[ip];
                    if (!p.valid) continue;
                    x.push_back(PHI_DEG[ip]);
                    y.push_back(p.bsa);
                    ey.push_back(std::max(1e-6, p.err));
                }
                if (!x.empty()) {
                    TGraphErrors* gr = new TGraphErrors((int)x.size(), x.data(), y.data(), nullptr, ey.data());
                    gr->SetMarkerStyle(20);
                    gr->SetMarkerSize(1.1);
                    gr->SetLineWidth(2);
                    gr->Draw("P SAME");
                }

                // Fitted curve (thin red line)
                if (cr.fit.status == 0 || cr.fit.ndf > 0) {
                    // draw via dense sampling to avoid any eval issues
                    const int NS=721;
                    std::vector<double> xd(NS), yd(NS);
                    for (int i=0;i<NS;++i){
                        double deg = double(i)*0.5;     // 0..360 in 0.5°
                        double rad = deg * (TWO_PI/360.0);
                        double denom = 1.0 + cr.fit.B1*std::cos(rad);
                        if (denom < EPS_DEN_EVAL) denom = EPS_DEN_EVAL;
                        double val = cr.fit.C + (cr.fit.A*std::sin(rad))/denom;
                        xd[i] = deg; yd[i] = val;
                    }
                    TGraph* gfit = new TGraph(NS, xd.data(), yd.data());
                    gfit->SetLineColor(kRed);
                    gfit->SetLineWidth(2);
                    gfit->Draw("L SAME");
                }

                // Panel annotation (Q² and -t)
                TLatex lab;
                lab.SetNDC(); lab.SetTextSize(0.040); lab.SetTextAlign(11);
                lab.SetTextFont(42);
                lab.DrawLatex(0.15, 0.94,
                    Form("Q^{2} #in [%.2g, %.2g],   -t #in [%.2g, %.2g]",
                         Q2_slice[ccol].first, Q2_slice[ccol].second,
                         t_slice[r].first,    t_slice[r].second));

                // Legend: keep well inside the pad
                TLegend* leg = new TLegend(0.50, 0.68, 0.90, 0.92);
                leg->SetBorderSize(1);
                leg->SetLineColor(kBlack);
                leg->SetFillColor(kWhite);
                leg->SetFillStyle(1001); 
                leg->SetTextFont(42);
                leg->SetTextSize(0.040);
                leg->AddEntry((TObject*)nullptr, Form("A = %.3f #pm %.3f",  cr.fit.A,  cr.fit.Aerr), "");
                leg->AddEntry((TObject*)nullptr, Form("B = %.3f #pm %.3f",  cr.fit.B1, cr.fit.B1err), "");
                leg->AddEntry((TObject*)nullptr, Form("C = %.3f #pm %.3f",  cr.fit.C,  cr.fit.Cerr), "");
                leg->Draw();
            }
        }

        std::ostringstream fout;
        fout << out_dir_plots << "/plot_bsa_" << period << "_xB_" << ix << ".png";
        c->SaveAs(fout.str().c_str());
        delete c;
    }
}

// ------------ helpers ------------
static inline bool inTenSix(const std::string& runTag) {
    return (runTag=="sp18_inb" || runTag=="sp18_out" ||
            runTag=="fa18_inb" || runTag=="fa18_inb_supp" || runTag=="fa18_out");
}

} // anon namespace

// =======================================================
// Public driver
// =======================================================
void compute_and_plot_bsa_helicity(
    const std::vector<std::string>& periods,
    const std::vector<std::string>& topologies,
    const std::vector<Binning>& binning_scheme,
    const std::map<std::string, TTree*>& dvcsDataTrees,
    const std::string& total_counts_json_path,
    const std::string& contamination_dir,
    const std::string& out_root_dir)
{
    namespace fs = std::filesystem;

    const auto xB_bins = uniqueRanges(binning_scheme, 'x');
    const auto Q2_bins = uniqueRanges(binning_scheme, 'Q');
    const auto t_bins  = uniqueRanges(binning_scheme, 't');

    GroupCounts groups;
    if (!load_total_counts(total_counts_json_path, groups)) {
        std::cerr<<"[bsa][ERROR] Failed to load total_counts json.\n"; return;
    }

    const auto PHI_RAD = phiCentersRad();

    const fs::path json_period_dir = fs::path(out_root_dir)/"jsons"/"BSA_fits";
    fs::create_directories(json_period_dir);

    std::map<std::string, std::map<std::tuple<int,int,int>, CellResult>> allPeriodCells;

    struct Acc { double Np=0, Nm=0; }; // effective counts sum
    std::map<std::tuple<int,int,int,int>, Acc> acc106; // (ix,iQ,it,ip)

    for (const auto& period : periods) {
        const std::string runTag = periodToRunTagKey(period);

        auto itG = groups.find(runTag);
        if (itG == groups.end()) {
            std::cerr<<"[bsa][WARN] total_counts has no group '"<<runTag<<"'. Skipping "<<period<<"\n";
            continue;
        }
        const auto& countsMap = itG->second;

        ContamMap contam;
        const fs::path contam_path = fs::path(contamination_dir)/("contamination_"+period+".json");
        if (!load_contam_for_period(contam_path.string(), contam)) {
            std::cerr<<"[bsa][WARN] No contamination file for "<<period<<". Assuming c=0.\n";
        }

        // polarization
        TTree* t = nullptr;
        auto itT = dvcsDataTrees.find(runTag);
        if (itT!=dvcsDataTrees.end()) t = itT->second;
        std::vector<PolStats> Pbin;
        PolStats Pavg{1.0,0.0,0};
        if (t){
            Pbin = compute_bin_polarization(t, xB_bins, Q2_bins, t_bins, topologies);
            double sumP=0; int nn=0;
            for (auto& ps : Pbin) if (ps.n>0) { sumP+=ps.P; nn++; }
            if (nn>0) { Pavg.P = sumP/double(nn); Pavg.n=nn; }
        }

        std::map<std::tuple<int,int,int>, CellResult> cells;

        // build cells
        for (int ix=0; ix<(int)xB_bins.size(); ++ix)
        for (int iQ=0; iQ<(int)Q2_bins.size(); ++iQ)
        for (int itb=0; itb<(int)t_bins.size(); ++itb) {
            CellResult result;
            result.points.resize(N_PHI_BINS);
            result.P_per_bin = true;
            result.P_used = (Pavg.P>0? Pavg.P : 1.0);

            for (int ip=0; ip<N_PHI_BINS; ++ip){
                const BinKey bk(ix,iQ,itb,ip);

                // raw counts
                auto itC = countsMap.find(bk);
                long long Np=0, Nm=0;
                if (itC!=countsMap.end()) { Np = itC->second.plus; Nm = itC->second.minus; }

                // contamination
                ContamBin cb;
                auto itCt = contam.find(bk);
                if (itCt!=contam.end()) cb = itCt->second;

                // corrected counts
                double Np_corr = Np*(1.0 - cb.c_plus);
                double Nm_corr = Nm*(1.0 - cb.c_minus);

                // polarization to use for this bin
                double P_here = result.P_used;
                if (!Pbin.empty()){
                    size_t idx = (((ix*(size_t)Q2_bins.size()+iQ)*(size_t)t_bins.size()+itb)*N_PHI_BINS + ip);
                    if (idx<Pbin.size() && Pbin[idx].n>0 && Pbin[idx].P>0.1) { P_here = Pbin[idx].P; }
                    else result.P_per_bin = false;
                }
                if (P_here<=0.0) { P_here = (Pavg.P>0? Pavg.P : 1.0); result.P_per_bin=false; }
                result.P_used = P_here;

                // scale by 1/P
                double Np_pol = Np_corr / P_here;
                double Nm_pol = Nm_corr / P_here;

                // --- Regularized BSA using Jeffreys prior (α = 0.5) ---
                const double alpha = 0.5;

                // Effective counts (pseudo-counts added to avoid ±1 with zero error)
                double Np_eff = std::max(0.0, Np_pol) + alpha;
                double Nm_eff = std::max(0.0, Nm_pol) + alpha;

                // Total and difference
                double S = Np_eff + Nm_eff;
                double D = Np_eff - Nm_eff;

                BSApt p; p.phi = PHI_RAD[ip];

                if (S > 0.0) {
                    p.bsa = D / S;

                    // Beta posterior variance mapped to A = 2p - 1
                    double a = Np_eff, b = Nm_eff;
                    double varA = 4.0 * (a * b) / ((a + b) * (a + b) * (a + b + 1.0));

                    // Safety floor to avoid absurdly small uncertainties
                    p.err = std::sqrt(std::max(varA, 1e-6));
                    p.valid = std::isfinite(p.bsa) && std::isfinite(p.err);
                } else {
                    p.bsa = 0.0;
                    p.err = 0.0;
                    p.valid = false;
                }

                result.points[ip] = p;

                // accumulate for combined 10.6 (using effective counts)
                if (inTenSix(runTag)) {
                    auto key4 = std::make_tuple(ix, iQ, itb, ip);
                    Acc& A = acc106[key4];
                    A.Np   += Np_eff;
                    A.Nm   += Nm_eff;
                }
            }

            result.fit = fit_cell(result.points);
            cells[std::make_tuple(ix,iQ,itb)] = std::move(result);
        }

        // write per-period JSON + plots
        const fs::path outP = json_period_dir/("BSA_fits_"+period+".json");
        write_period_bsa_json(outP.string(), N_PHI_BINS, xB_bins, Q2_bins, t_bins, cells);

        const fs::path plots_dir = fs::path(out_root_dir)/"bsa_plots"/periodToRunTagKey(period);
        std::error_code ec; fs::create_directories(plots_dir, ec);
        plot_cells_for_period(period, binning_scheme, xB_bins, Q2_bins, t_bins, cells, plots_dir.string());

        allPeriodCells[period] = std::move(cells);
    }

    // ---- write “all periods” file ----
    write_all_periods_json((fs::path(out_root_dir)/"jsons"/"BSA_fits_all_periods.json").string(),
                           allPeriodCells, N_PHI_BINS, xB_bins, Q2_bins, t_bins);

    // ---- 10.6 GeV statistical combination ----
    std::map<std::tuple<int,int,int>, CellResult> combCells;
    for (int ix=0; ix<(int)xB_bins.size(); ++ix)
    for (int iQ=0; iQ<(int)Q2_bins.size(); ++iQ)
    for (int itb=0; itb<(int)t_bins.size(); ++itb) {
        CellResult cr; cr.points.resize(N_PHI_BINS);
        cr.P_used = 1.0; cr.P_per_bin = false;
        for (int ip=0; ip<N_PHI_BINS; ++ip){
            auto itA = acc106.find(std::make_tuple(ix,iQ,itb,ip));
            BSApt p; p.phi = phiCentersRad()[ip]; p.valid=false;
            if (itA!=acc106.end()){
                const auto& A = itA->second;

                const double alpha = 0.5;
                double Np_eff = std::max(0.0, A.Np); // already had alpha in accumulation
                double Nm_eff = std::max(0.0, A.Nm);

                double S = Np_eff + Nm_eff;
                double D = Np_eff - Nm_eff;

                if (S > 0.0) {
                    p.bsa = D / S;

                    // Conservative uncertainty (as if a,b already included α from each period)
                    double a = Np_eff, b = Nm_eff;
                    double varA = 4.0 * (a * b) / ((a + b) * (a + b) * (a + b + 1.0));
                    p.err = std::sqrt(std::max(varA, 1e-6));
                    p.valid = std::isfinite(p.bsa) && std::isfinite(p.err);
                }
            }
            cr.points[ip]=p;
        }
        cr.fit = fit_cell(cr.points);
        combCells[std::make_tuple(ix,iQ,itb)] = std::move(cr);
    }

    // write combined JSON (uses 10.6 with a dot)
    {
        std::ofstream ofs((fs::path(out_root_dir)/"jsons"/"BSA_fits_combined_10.6.json").string());
        if (ofs){
            ofs<<std::fixed<<std::setprecision(8);
            ofs<<"{\n";
            ofs<<"  \"combined\": true,\n";
            ofs<<"  \"periods_used\": [\"DVCS_Sp18_inb\",\"DVCS_Sp18_out\",\"DVCS_Fa18_inb_supp\",\"DVCS_Fa18_inb\",\"DVCS_Fa18_out\"],\n";
            ofs<<"  \"binning_meta\": {\"phi_bins\": "<<N_PHI_BINS<<", \"xB_bins\": "<<xB_bins.size()<<", \"Q2_bins\": "<<Q2_bins.size()<<", \"t_bins\": "<<t_bins.size()<<"},\n";
            ofs<<"  \"bins\": {\n";
            bool first=true;
            for (const auto& kv : combCells){
                if (!first) ofs<<",\n"; first=false;
                int ix,iQ,it; std::tie(ix,iQ,it)=kv.first;
                const auto& cr = kv.second;
                ofs<<"    \"("<<ix<<","<<iQ<<","<<it<<")\": {";
                ofs<<"\"phi\":["; for (size_t i=0;i<cr.points.size();++i){ if(i)ofs<<","; ofs<<cr.points[i].phi; } ofs<<"],";
                ofs<<"\"bsa\":["; for (size_t i=0;i<cr.points.size();++i){ if(i)ofs<<","; ofs<<cr.points[i].bsa; } ofs<<"],";
                ofs<<"\"bsa_err\":["; for (size_t i=0;i<cr.points.size();++i){ if(i)ofs<<","; ofs<<cr.points[i].err; } ofs<<"],";
                ofs<<"\"fit\":{"
                      "\"A\":{\"val\":"<<cr.fit.A<<",\"err\":"<<cr.fit.Aerr<<"},"
                      "\"B1\":{\"val\":"<<cr.fit.B1<<",\"err\":"<<cr.fit.B1err<<"},"
                      "\"B2\":{\"val\":0.0,\"err\":0.0},"
                      "\"C\":{\"val\":"<<cr.fit.C<<",\"err\":"<<cr.fit.Cerr<<"},"
                      "\"chi2\":"<<cr.fit.chi2<<",\"ndf\":"<<cr.fit.ndf<<",\"status\":"<<cr.fit.status<<"}";
                ofs<<"}";
            }
            ofs<<"\n  }\n}\n";
        }
    }

    // plots for combined 10.6 (directory name matches the dot style)
    const fs::path plots_comb106 = fs::path(out_root_dir)/"bsa_plots"/"10.6_combined";
    std::error_code ec; fs::create_directories(plots_comb106, ec);
    plot_cells_for_period("RGA_10.6_combined", binning_scheme, xB_bins, Q2_bins, t_bins, combCells, plots_comb106.string());
}