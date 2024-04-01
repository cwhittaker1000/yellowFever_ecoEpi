/* IBM Unica Page Tagging Script v2.2.0
 *
 * Licensed Materials - Property of IBM (c) Copyright IBM Corporation 2004,2011.
 * U.S. Government Users Restricted Rights: Use, duplication or disclosure restricted by GSA ADP Schedule Contract with IBM Corp.
 */

var NTPT_IMGSRC = 'http://pftag.scholarlyiq.com/siqpagetag.gif';

var NTPT_FLDS = {};
NTPT_FLDS.lc = true; // Document location
NTPT_FLDS.rf = true; // Document referrer
NTPT_FLDS.rs = true; // User's screen resolution
NTPT_FLDS.cd = true; // User's color depth
NTPT_FLDS.ln = true; // Browser language
NTPT_FLDS.tz = true; // User's timezone
NTPT_FLDS.jv = true; // Browser's Java support
NTPT_FLDS.ck = true; // Cookies
NTPT_FLDS.iv = false; // Initial view

var NTPT_MAXTAGWAIT = 1.0; // Max delay (secs) on link-tags and submit-tags

// Optional variables:
var NTPT_HTTPSIMGSRC = 'https://pftag.scholarlyiq.com/siqpagetag.gif';
var NTPT_GLBLEXTRA = '';
var NTPT_GLBLREFTOP = false;
var NTPT_GLBLCOOKIES = [ ];

/*** END OF USER-CONFIGURABLE VARIABLES ***/
(function(){var a=document,r=window;function Q(value){return((typeof(value)==="string")&&(value!==""));};function x5(k){return(encodeURIComponent?encodeURIComponent(k):escape(k));};function mr(k){return(decodeURIComponent?decodeURIComponent(k):unescape(k));};function p(B,k,wF,ip){var H="",T;H=B+'='+x5(k)+";";if(ip){H+=" domain="+ip+";";}if(wF>0){T=new Date();T.setTime(T.getTime()+(wF*1000));H+=" expires="+T.toGMTString()+";";}H+=" path=/";a.cookie=H;};function W(B){var U,M,b,H;if(Q(B)){U=B+"=";H=a.cookie;if(H.length>0){b=ZO(U,H,0);if(b!==-1){b+=U.length;M=H.indexOf(";",b);if(M===-1){M=H.length;}return mr(H.substring(b,M));};};}return null;};function ZO(U,H,SL){var b=H.indexOf(U,SL);if(b<0){return-1;}else if((b===0)||((b>1)&&(H.substring(b-2,b)==="; "))){return b;}return ZO(U,H,(b+1));};function or(wS){var L="",G;for(G in wS){if(Q(wS[G])){if(L!==""){L+=";";}L+=G+"="+wS[G];};}return L;};function d1(f){var nO='ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789';if(f<62){return nO.charAt(f);}return(d1(Math.floor(f/62))+nO.charAt(f%62));};function Sl(){var JP="",CB=new Date(),Sx;for(Sx=0;Sx<11;Sx+=1){JP+=d1(Math.round(Math.random()*61));}return(JP+"-"+d1(CB.getTime()));};function BR(m,nT){return(m+(((m==='')||((nT==='')||(nT.substring(0,1)==='&')))?'':'&')+nT);};function jJ(){var cQ=new Date();return(cQ.getTime()+'.'+Math.floor(Math.random()*1000));};function h(oa,KN){return(typeof(r[oa])!=="undefined")?r[oa]:KN;};(function(){var sZ=h('NTPT_GLBLCOOKIES',null),kU=h('NTPT_PGCOOKIES',null),ED=h('NTPT_SET_IDCOOKIE',true),K=h('NTPT_IDCOOKIE_NAME',"SaneID"),hE=h('NTPT_IDCOOKIE_DOMAIN',null),va=h('NTPT_IDCOOKIE_EXPIRE',155520000),yR=h('NTPT_SET_SESSION_COOKIE',false),PN=h('NTPT_SESSION_COOKIE_NAME',"NetInsightSessionID"),vL=h('NTPT_SESSION_COOKIE_DOMAIN',null),tN=h('NTPT_HTTPSIMGSRC',''),Nt=h('NTPT_PGREFTOP',h('NTPT_GLBLREFTOP',false)),_Q=h('NTPT_NOINITIALTAG',false),NY=h('NTPT_MAXTAGWAIT',1.0),E2=NTPT_IMGSRC,c=NTPT_FLDS,hP=null,i=null,q=null,D=null,A=[],d={},C=new Array(10),z;for(z=0;z<C.length;z+=1){C[z]=null;}function v(B,JS){if(typeof(JS)!=="undefined"){A[B]=JS.toString();}};function T9(B){A[B]='';};function q1(S){var sa,I,aQ;if(Q(S)){sa=S.split('&');for(aQ=0;aQ<sa.length;aQ+=1){I=sa[aQ].split('=');if(I.length===2){v(I[0],mr(I[1]));};};}};function ZH(S){var Os='',G;q1(h('NTPT_GLBLEXTRA',''));if(!q){q1(h('NTPT_PGEXTRA',''));}q1(S);for(G in A){if(Q(A[G])){Os=BR(Os,G+'='+x5(A[G]));};}return Os;};function kp(){var G;d.A=[];for(G in A){if(Q(A[G])){d.A[G]=A[G];};}};function JF(){var G;A=[];for(G in d.A){if(Q(d.A[G])){A[G]=d.A[G];};}};function A4(P,e,w){if(C[P]!==null){C[P].onload=e;C[P].onerror=e;C[P].onabort=e;}return setTimeout(e,(w*1000));};function pn(N,R){if(!Q(N)){return;}z=((z+1)%C.length);if(C[z]===null){C[z]=new Image(1,1);}C[z].src=N+'?'+R;};function jy(S){var N,R;if((tN!=='')&&(a.location.protocol==='https:')){N=tN;}else{N=E2;}R=ZH(S);pn(N,R);JF();};function Zt(S){v('ets',jJ());jy(S);return true;};function XT(Z){if(hP){clearTimeout(hP);hP=null;}if(q){var _=q;q=null;if(Z){_.click();_.onclick=_.TD;}else{r.location.href=_.href;};}};function lD(t,S,w){var o,Z;if(!t||!t.href){return true;}if(q){return false;}q=t;if(c.lc){v('lc',t.href);}if(c.rf){if(!Nt||!top||!top.document){v('rf',a.location);};}Zt(S);if(w){o=w;}else{o=NY;}if(o>0&&(q.target==""||q.target==r.name)){if(q.click){q.TD=q.onclick;q.onclick=null;Z=true;}else{Z=false;}hP=A4(z,function(){XT(Z);},o);return false;}q=null;return true;};function UV(){if(i){clearTimeout(i);i=null;}if(D){var Og=D;D=null;Og.submit();Og.onsubmit=Og.Pe;}};function Fs(Y,S,w){var o;if(!Y||!Y.submit){return true;}if(D){return false;}D=Y;Zt(S);if(w){o=w;}else{o=NY;}if(o>0){Y.Pe=Y.onsubmit;Y.onsubmit=null;i=A4(z,function(){UV();},o);return false;}D=null;return true;};function iS(){var et;if(Nt&&top&&top.document){et=top.document.referrer;}else{et=a.referrer;}v('rf',et);};function j9(){var O;if(navigator.language){O=navigator.language;}else if(navigator.userLanguage){O=navigator.userLanguage;}else{O='';}if(O.length>2){O=O.substring(0,2);}O=O.toLowerCase();v('ln',O);};function PQ(){var J,cQ=new Date(),l=cQ.getTimezoneOffset(),V;J='GMT';if(l!==0){if(l>0){J+=' -';}else{J+=' +';}l=Math.abs(l);V=Math.floor(l/60);l-=V*60;if(V<10){J+='0';}J+=V+':';if(l<10){J+='0';}J+=l;}v('tz',J);};function kt(){var X=[],NB=false,Pq='ck',s,H,g,G;v('js','1');v('ts',jJ());if(c.lc){v('lc',a.location);}if(c.rf){iS();}if(self.screen){if(c.rs){v('rs',self.screen.width+'x'+self.screen.height);}if(c.cd){v('cd',self.screen.colorDepth);};}if(c.ln){j9();}if(c.tz){PQ();}if(c.jv){v('jv',navigator.javaEnabled()?'1':'0');}if(yR&&!W(PN)){p(PN,1,0,vL);if(c.iv&&W(PN)){v('iv','1');};}if(c.ck){if(sZ){for(s=0;s<sZ.length;s+=1){X[sZ[s]]="";};}if(kU){for(s=0;s<kU.length;s+=1){X[kU[s]]="";};}for(G in X){if(typeof(X[G])==="string"){H=W(G);if(H){X[G]=H;};};}if(ED){H=W(K);if(H){X[K]=H;NB=true;};}g=or(X);if(g!==""){v(Pq,g);};}kp();if(!_Q){jy('');}T9('iv');kp();if(ED&&!NB){H=W(K);if(!H){H=Sl();p(K,H,va,hE);if(c.ck&&W(K)){X[K]=H;g=or(X);if(g!==""){v(Pq,g);kp();};};};}};function bx(B,k){var m,Mz,s,Jw;m=a.location.search.substring(1);Jw=B+"="+k;Mz=m.split("&");for(s=0;s<Mz.length;s+=1){if(Mz[s]===Jw){return true;};}return false;};function UX(){var U0=h("NTPT_EM_COOKIE_PREFIX","NetInsightEM"),nS=U0+"Session",RH=U0,sx=h("NTPT_EM_COOKIE_TIMEOUT",1800),RS="emsgpcat",MO="1",kz="1";if(W(nS)||W(RH)||bx(RS,MO)){p(nS,kz,0,hE);p(RH,kz,sx,hE);v(RS,MO);return true;}return false;};function Y1(){if(c.gsme){return UX();}return true;};if(Y1()){r.ntptAddPair=function(F,j){return v(F,j);};r.ntptDropPair=function(F){return T9(F);};r.ntptEventTag=function(F){return Zt(F);};r.ntptLinkTag=function(F,j,nY){return lD(F,j,nY);};r.ntptSubmitTag=function(F,j,nY){return Fs(F,j,nY);};kt();}else{r.ntptAddPair=r.ntptDropPair=r.ntptEventTag=r.ntptLinkTag=r.ntptSubmitTag=function(){return true;};}}());}());
