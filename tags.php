<?
 function search($text_search){

  // Tidy inputs
  $text_search = explode(" ", $text_search);
  $strip = array(",", " ", "#");
  $i = 0;
  while($i<count($text_search))
  {
    $text_search[$i] = str_replace($strip,"",$text_search[$i]);
    $i++;
  }

  // Load data
  $buffer = file('tags.list');

  // Search
  $i = 0;
  while($i<count($text_search)) {
   $j = 0;
   $buffer2 = array();
   while($j<count($buffer)) {
      if ( strripos($buffer[$j], $text_search[$i]) !== false ) {
         array_push($buffer2, $buffer[$j]);
      }
      $j++;
   }
   $buffer = $buffer2;
   $i++;
  }

  // Tidy Output
  $j = 0;
  $found = array();
  while($j<count($buffer)) {
    $buffer[$j] = explode(" ", $buffer[$j]);
    array_push($found, $buffer[$j][0]);
    $j++;
  }

  return $found;
 }
?>
<head>
<title>Tag search</title>
<meta http-equiv="content-type" content="text/html; charset=utf-8" />
<link rel="shortcut icon" href="/chebfun/images/chebicon_32x32.ico" type="image/x-icon" />
<meta name="keywords" content="Chebfun, object-oriented matlab, functions" />
<link href="/chebfun/css/style.css" rel="stylesheet" type="text/css" />
<link href="/chebfun/css/mainmenu.css" rel="stylesheet" type="text/css" />
<link rel="image_src" href="/chebfun/images/cheblogo_c_small.png"/>
<link href="/chebfun/css/print.css" rel="stylesheet" type="text/css" media="print" />
<script type="text/javascript" src="https://apis.google.com/js/plusone.js">
  {lang: 'en-GB'}
</script>
<script type="text/javascript"
  src="https://c328740.ssl.cf1.rackcdn.com/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML">
</script>
</head>

<body>
<div id="headimg"></div>
<div id="wrapper"> <!--WRAPPER-->	<div id="myheader"> <!-- LOGO & HEADER -->
		<a href="/chebfun/"><img src="/chebfun/images/cheblogow_small.png" alt="Chebfun Logo" class="logoimg" width="300px" /></a>
	  <div style="float:right; width:300px; margin-top:5px;">


        <div style="float:right; width:300px; margin-top:1px; margin-bottom:2px">
    		<a class="headerA" style="margin-top:0px;" href="http://www.ox.ac.uk" title="University of Oxford home page">Oxford&nbsp;University</a><br/>
    		<a class="headerB" href="/" title="Mathematical Institute home page">Mathematical&nbsp;Institute</a>       
	    </div>


        <!-- Place this tag where you want the +1 button to render -->
        <div style="float:right; width:63px; margin-left:2px; margin-right:5px"><g:plusone size="medium"></g:plusone></div>

	    <div style="width:83px; height:20px; overflow:none; float:right; margin-left:5px;" >
		<iframe src="http://www.facebook.com/plugins/like.php?href=http://www.chebfun.org&layout=button_count&show_faces=false&width=80&action=like&colorscheme=light" 
        	scrolling="no" frameborder="0"	style="border:none; width:95px; height:20px"></iframe>
	    </div> 
<!--            <div id="fbsharediv" style="width:65px; height=25px; overflow=false; padding-top:1px;">
                <a name="fb_share" type="button" share_url="http://www.maths.ox.ac.uk/chebfun/" href="http://www.facebook.com/sharer.php">Share</a><script src="http://static.ak.fbcdn.net/connect.php/js/FB.Share" type="text/javascript"></script>
            </div>-->
<!--        <div id="tweetdiv" style="margin-left:7px;"> 
		<a href="http://twitter.com/share" class="twitter-share-button" data-url="http://www.chebfun.org" data-text="#Chebfun is an open-source software system for #MATLAB which extends methods of computation for numbers to functions." data-count="none">Tweet</a><script type="text/javascript" src="http://platform.twitter.com/widgets.js"></script>
        </div>-->

        <div id="tweetdiv" style="margin-left:7px; width:61px;"> 
            <a href='https://twitter.com/share?original_referer=http
                <img alt='' border='0' src='/chebfun/images/tweet2.png' name='tweet' style="width:61px; height:20px;"/>
            </a>
        </div>

        <div style="width:61px; float:right;">
            <a href='http://twitter.com/chebfun' onmouseover="document.follow.src='/chebfun/images/follow.png'" onmouseout="document.follow.src='/chebfun/images/follow2.png'">
                <img alt='' border='0' src='/chebfun/images/follow2.png' name='follow' style="width:61px; height:20px;"/>
            </a>
        </div>

      </div>
	</div> <!-- END LOGO & HEADER -->


 
<nav> <!-- END MAIN MENU -->
<div id="main-menu">
    <b class="rtop"><b class="r1"></b><b class="r2"></b><b class="r3"></b><b class="r4"></b></b>
	<ul>
		<li><a href="/chebfun/">home</a></li>
   		<li><a href="/chebfun/about/">about</a>
            <ul>
    	        <li><a href="/chebfun/about/">The software</a></li>
            	<li><a href="/chebfun/about/">The team</a></li>
            	<li><a href="/chebfun/and_beyond/">2012 workshop</a></li>
        	</ul>
        </li>
<!--			<li><a href="/chebfun/news/">news</a></li> -->
		<li><a href="/chebfun/download/">download</a>
		    <ul>
		    	<li><a href="/chebfun/download/">stable release</a></li>
		    	<li><a href="/chebfun/nightly/">nightly release</a></li>
    		</ul>
        </li>
        <li><a href="/chebfun/guide/">documentation</a>
            <ul>
                <li><a href="/chebfun/guide/">Chebfun guide</a></li>
                <li><a href="/chebfun/chebfun2/guide/">Chebfun2 guide</a></li>
            </ul>
        </li>
	    <li><a href="/chebfun/examples/">examples</a>
            <ul>
		    	<li><a href="/chebfun/examples/">Chebfun examples</a></li>
		    	<li><a href="/chebfun/examples/tags.php?query=chebfun2">Chebfun2 examples</a></li>
		    	<li><a href="/chebfun/examples/new/">new examples</a></li>
		    	<li><a href="/chebfun/examples/submit.shtml">submit an example</a></li>
    		</ul>
        </li>
		<li><a href="/chebfun/faq/">faq</a></li>
        <li><a href="/chebfun/chebfun2/">Chebfun2</a>
            <ul>
                <li><a href="/chebfun/chebfun2/">about</a></li>
<!--                <li><a href="/chebfun/chebfun2/">download</a></li>-->
                <li><a href="/chebfun/chebfun2/guide/">guide</a></li>
                <li><a href="/chebfun/examples/tags.php?query=chebfun2">examples</a></li>
            </ul>
        </li>
		<li><a href="/chebfun/publications/">publications</a></li> 
<!--			<li><a href="https://svn.maths.ox.ac.uk/trac/chebfun/wiki/ChebfunPublications">publications</a></li> -->
  	    <li><a href="https://svn.maths.ox.ac.uk/trac/chebfun/">bugs</a>
		    <ul>
		       	<li><a href="https://svn.maths.ox.ac.uk/trac/chebfun/newticket">submit a bug report</a></li>
    		</ul>
        </li>
        <li><a href="https://svn.maths.ox.ac.uk/trac/chebfun/wiki/DeveloperZone">develop</a>
			<ul>
				<li><a href="https://svn.maths.ox.ac.uk/trac/chebfun/browser/trunk">svn repository</a></li>
				<li><a href="https://svn.maths.ox.ac.uk/trac/chebfun/wiki/">Trac</a></li>
			</ul>
	   </li>
<!--		<li><a href="/chebfun/contact/">contact</a>
			<ul>
				<lil><a href="/chebfun/lists/">mailing lists</a></li>
			</ul>
		</li>
-->
    	<li><a href="/chebfun/search/">search</a></li>
	</ul>
	<b class="rbottom"><b class="r4"></b><b class="r3"></b><b class="r2"></b><b class="r1"></b></b>
</div> 
</nav> <!-- END MAIN MENU -->

<h2 class="bigger">Tag search</h2>
<div id="mycontent"> <!--CONTENT-->
<form action="tags.php" method="get">
search: <input type="text" name="query" />
<input type="submit" />
</form><?php
  $string = $_GET["query"]; 

  $found = search($string);
  $i = 0; 
  if (count($string) !== 0) {
    if (count($found) == 0) {
      print("<h2 class='bigger' style='margin-left:-15pt; margin-bottom:-10pt;'>No search results for '$string'.</h2>");
    }
    else {
      print("<h2 class='bigger' style='margin-left:-15pt; margin-bottom:-10pt;'>Search results for '$string':</h2>");
    }
  }
  while($i<count($found))
  {
    $found[$i] = explode("/", $found[$i]);
    $dir = $found[$i][0];
    $file = $found[$i][1];
    $url = "http://www2.maths.ox.ac.uk/chebfun/examples/$dir/html/$file.shtml";
    $mfile = @fopen("http://www2.maths.ox.ac.uk/chebfun/examples/$dir/$file.m", "r");
    $title = substr(fgets($mfile, 4096), 3);
    fclose($mfile);
    print("<br/><a href='$url'>$title</a> : ($dir/$file)");
    $i++;
  }
?>
<h2 class="bigger" style="margin-left:-15pt; margin-top:10pt; margin-bottom:-10pt;">Tag cloud</h2><br/>

<span style="width:400px"><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=Xmas">Xmas</a></b><b style="font-size:100%"> </b><b style="font-size:104%"><a href="/chebfun/examples/tags.php?query=spline">spline</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=qrdecomposition">qrdecomposition</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=centrallimittheorem">centrallimittheorem</a></b><b style="font-size:100%"> </b><b style="font-size:104%"><a href="/chebfun/examples/tags.php?query=Chebyshevpolynomial">Chebyshevpolynomial</a></b><b style="font-size:100%"> </b><b style="font-size:136%"><a href="/chebfun/examples/tags.php?query=jumpconditions">jumpconditions</a></b><b style="font-size:100%"> </b><b style="font-size:104%"><a href="/chebfun/examples/tags.php?query=impulse">impulse</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=concatentation">concatentation</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=curve intersections">curve intersections</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=resonance">resonance</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=splitting">splitting</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=continuation">continuation</a></b><b style="font-size:100%"> </b><b style="font-size:215%"><a href="/chebfun/examples/tags.php?query=fun">fun</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=pretty">pretty</a></b><b style="font-size:100%"> </b><b style="font-size:460%"><a href="/chebfun/examples/tags.php?query=linearODE">linearODE</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=Resultant">Resultant</a></b><b style="font-size:100%"> </b><b style="font-size:136%"><a href="/chebfun/examples/tags.php?query=norm">norm</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=waves">waves</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=resolution">resolution</a></b><b style="font-size:100%"> </b><b style="font-size:238%"><a href="/chebfun/examples/tags.php?query=geometry">geometry</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=polygon">polygon</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=sinewaves">sinewaves</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=fast transform">fast transform</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=binomial">binomial</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=signaprocessing">signaprocessing</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=quadrature_weights">quadrature_weights</a></b><b style="font-size:100%"> </b><b style="font-size:340%"><a href="/chebfun/examples/tags.php?query=nonlinearODE">nonlinearODE</a></b><b style="font-size:100%"> </b><b style="font-size:281%"><a href="/chebfun/examples/tags.php?query=rational">rational</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=parameterized">parameterized</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=fieldofvalues">fieldofvalues</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=soccer">soccer</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=love">love</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=chaos">chaos</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=contourintegral">contourintegral</a></b><b style="font-size:100%"> </b><b style="font-size:104%"><a href="/chebfun/examples/tags.php?query=dynamical systems">dynamical systems</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=implicitfunction">implicitfunction</a></b><b style="font-size:100%"> </b><b style="font-size:340%"><a href="/chebfun/examples/tags.php?query=piecewise">piecewise</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=Pachon">Pachon</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=fouriertransform">fouriertransform</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=corners">corners</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=Reimann">Reimann</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=Bouncing">Bouncing</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=Leja">Leja</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=fractal">fractal</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=ode">ode</a></b><b style="font-size:100%"> </b><b style="font-size:136%"><a href="/chebfun/examples/tags.php?query=linearPDE">linearPDE</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=highorder">highorder</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=nonnormality">nonnormality</a></b><b style="font-size:100%"> </b><b style="font-size:104%"><a href="/chebfun/examples/tags.php?query=phase portrait">phase portrait</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=curve">curve</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=contourc">contourc</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=quadrature_nodes">quadrature_nodes</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=Divergentseries">Divergentseries</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=Lorenz">Lorenz</a></b><b style="font-size:100%"> </b><b style="font-size:136%"><a href="/chebfun/examples/tags.php?query=convolution">convolution</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=exponents">exponents</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=erosion">erosion</a></b><b style="font-size:100%"> </b><b style="font-size:281%"><a href="/chebfun/examples/tags.php?query=optimization">optimization</a></b><b style="font-size:100%"> </b><b style="font-size:104%"><a href="/chebfun/examples/tags.php?query=coefficients">coefficients</a></b><b style="font-size:100%"> </b><b style="font-size:136%"><a href="/chebfun/examples/tags.php?query=pole">pole</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=Green's theorem">Green's theorem</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=arclength">arclength</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=polynomials">polynomials</a></b><b style="font-size:100%"> </b><b style="font-size:104%"><a href="/chebfun/examples/tags.php?query=best">best</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=Dawson">Dawson</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=entire">entire</a></b><b style="font-size:100%"> </b><b style="font-size:164%"><a href="/chebfun/examples/tags.php?query=deltafunction">deltafunction</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=metastability">metastability</a></b><b style="font-size:100%"> </b><b style="font-size:136%"><a href="/chebfun/examples/tags.php?query=trajectory">trajectory</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=specialfunctions">specialfunctions</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=integration">integration</a></b><b style="font-size:100%"> </b><b style="font-size:104%"><a href="/chebfun/examples/tags.php?query=residue">residue</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=resolvent">resolvent</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=physics">physics</a></b><b style="font-size:100%"> </b><b style="font-size:104%"><a href="/chebfun/examples/tags.php?query=integral">integral</a></b><b style="font-size:100%"> </b><b style="font-size:136%"><a href="/chebfun/examples/tags.php?query=quantum">quantum</a></b><b style="font-size:100%"> </b><b style="font-size:136%"><a href="/chebfun/examples/tags.php?query=ODEsystem">ODEsystem</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=Hermite">Hermite</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=localextrema">localextrema</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=Clenshaw_Curtis_quadrature">Clenshaw_Curtis_quadrature</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=reactiondiffusion">reactiondiffusion</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=integrodifferential">integrodifferential</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=normaldistribution">normaldistribution</a></b><b style="font-size:100%"> </b><b style="font-size:104%"><a href="/chebfun/examples/tags.php?query=wikipedia">wikipedia</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=special function">special function</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=Bezout">Bezout</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=line integral">line integral</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=CarrierPearson">CarrierPearson</a></b><b style="font-size:100%"> </b><b style="font-size:136%"><a href="/chebfun/examples/tags.php?query=nonstandard">nonstandard</a></b><b style="font-size:100%"> </b><b style="font-size:104%"><a href="/chebfun/examples/tags.php?query=conformalmap">conformalmap</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=GramSchmidt">GramSchmidt</a></b><b style="font-size:100%"> </b><b style="font-size:104%"><a href="/chebfun/examples/tags.php?query=FFT">FFT</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=Bernoulli">Bernoulli</a></b><b style="font-size:100%"> </b><b style="font-size:492%"><a href="/chebfun/examples/tags.php?query=Chebfun2">Chebfun2</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=volterra">volterra</a></b><b style="font-size:100%"> </b><b style="font-size:104%"><a href="/chebfun/examples/tags.php?query=Gauss_quadrature">Gauss_quadrature</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=hyperfunction">hyperfunction</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=Maxwell">Maxwell</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=minsamples">minsamples</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=exponentialdistribution">exponentialdistribution</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=Duffing oscillator">Duffing oscillator</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=Hilbert">Hilbert</a></b><b style="font-size:100%"> </b><b style="font-size:238%"><a href="/chebfun/examples/tags.php?query=probability">probability</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=Chebyshevpoints">Chebyshevpoints</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=dual">dual</a></b><b style="font-size:100%"> </b><b style="font-size:136%"><a href="/chebfun/examples/tags.php?query=linear">linear</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=random">random</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=mode">mode</a></b><b style="font-size:100%"> </b><b style="font-size:215%"><a href="/chebfun/examples/tags.php?query=interpolation">interpolation</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=parallelogram law">parallelogram law</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=angle">angle</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=chebsnake">chebsnake</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=Ball">Ball</a></b><b style="font-size:100%"> </b><b style="font-size:104%"><a href="/chebfun/examples/tags.php?query=contour">contour</a></b><b style="font-size:100%"> </b><b style="font-size:340%"><a href="/chebfun/examples/tags.php?query=eigenvalues">eigenvalues</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=Rank">Rank</a></b><b style="font-size:100%"> </b><b style="font-size:164%"><a href="/chebfun/examples/tags.php?query=astrophysics">astrophysics</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=discontinuity">discontinuity</a></b><b style="font-size:100%"> </b><b style="font-size:190%"><a href="/chebfun/examples/tags.php?query=rootfinding">rootfinding</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=heaviside">heaviside</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=Lanczos">Lanczos</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=Beta distribution">Beta distribution</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=normal distribution">normal distribution</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=Allen-Cahn">Allen-Cahn</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=stabilityregions">stabilityregions</a></b><b style="font-size:100%"> </b><b style="font-size:104%"><a href="/chebfun/examples/tags.php?query=analytic continuation">analytic continuation</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=AllenCahn">AllenCahn</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=snowflake>">snowflake></a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=divgradcurl">divgradcurl</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=deltafunctions">deltafunctions</a></b><b style="font-size:100%"> </b><b style="font-size:104%"><a href="/chebfun/examples/tags.php?query=mean">mean</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=compact">compact</a></b><b style="font-size:100%"> </b><b style="font-size:104%"><a href="/chebfun/examples/tags.php?query=variance">variance</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=l2approximation">l2approximation</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=periodic">periodic</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=exoticBCs">exoticBCs</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=fredholm">fredholm</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=automaticdifferention">automaticdifferention</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=camorbit">camorbit</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=bvp">bvp</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=Snell">Snell</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=basis">basis</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=Bspline">Bspline</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=Legendre">Legendre</a></b><b style="font-size:100%"> </b><b style="font-size:136%"><a href="/chebfun/examples/tags.php?query=phase portraits">phase portraits</a></b><b style="font-size:100%"> </b><b style="font-size:104%"><a href="/chebfun/examples/tags.php?query=stability">stability</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=affinetransformations">affinetransformations</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=conditionnumber">conditionnumber</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=fractional">fractional</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=birthday">birthday</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=100-Digit Challenge">100-Digit Challenge</a></b><b style="font-size:100%"> </b><b style="font-size:376%"><a href="/chebfun/examples/tags.php?query=complex">complex</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=shapeanaylsis">shapeanaylsis</a></b><b style="font-size:100%"> </b><b style="font-size:164%"><a href="/chebfun/examples/tags.php?query=delta">delta</a></b><b style="font-size:100%"> </b><b style="font-size:136%"><a href="/chebfun/examples/tags.php?query=exponential">exponential</a></b><b style="font-size:100%"> </b><b style="font-size:104%"><a href="/chebfun/examples/tags.php?query=intersection">intersection</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=bivariate polynomials">bivariate polynomials</a></b><b style="font-size:100%"> </b><b style="font-size:260%"><a href="/chebfun/examples/tags.php?query=linearalgebra">linearalgebra</a></b><b style="font-size:100%"> </b><b style="font-size:104%"><a href="/chebfun/examples/tags.php?query=bestapproximation">bestapproximation</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=squareroot">squareroot</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=Higham">Higham</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=Bernstein">Bernstein</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=pade">pade</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=Pade">Pade</a></b><b style="font-size:100%"> </b><b style="font-size:104%"><a href="/chebfun/examples/tags.php?query=Chebyshev">Chebyshev</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=Caratheodory">Caratheodory</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=transient">transient</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=determinant">determinant</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=Chebyshevcoefficients">Chebyshevcoefficients</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=median">median</a></b><b style="font-size:100%"> </b><b style="font-size:104%"><a href="/chebfun/examples/tags.php?query=surfaces">surfaces</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=gamma">gamma</a></b><b style="font-size:100%"> </b><b style="font-size:104%"><a href="/chebfun/examples/tags.php?query=subdivision">subdivision</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=blackscholes">blackscholes</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=perimeter">perimeter</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=Cauchy">Cauchy</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=gaussian">gaussian</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=backslash">backslash</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=Bessel">Bessel</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=parametric">parametric</a></b><b style="font-size:100%"> </b><b style="font-size:104%"><a href="/chebfun/examples/tags.php?query=spike">spike</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=Maxwell distribution">Maxwell distribution</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=Gibbsphenomenon">Gibbsphenomenon</a></b><b style="font-size:100%"> </b><b style="font-size:376%"><a href="/chebfun/examples/tags.php?query=roots">roots</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=mulivriable calculus">mulivriable calculus</a></b><b style="font-size:100%"> </b><b style="font-size:190%"><a href="/chebfun/examples/tags.php?query=IVP">IVP</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=bump">bump</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=FundamentalTheoremOfAlgebra">FundamentalTheoremOfAlgebra</a></b><b style="font-size:100%"> </b><b style="font-size:104%"><a href="/chebfun/examples/tags.php?query=2D">2D</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=singular_integral">singular_integral</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=whitecurves">whitecurves</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=system">system</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=heatequation">heatequation</a></b><b style="font-size:100%"> </b><b style="font-size:104%"><a href="/chebfun/examples/tags.php?query=approximation">approximation</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=Gibbs">Gibbs</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=Snell's Law">Snell's Law</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=multiple_integrals">multiple_integrals</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=Domain">Domain</a></b><b style="font-size:100%"> </b><b style="font-size:104%"><a href="/chebfun/examples/tags.php?query=Newton">Newton</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=sound">sound</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=koch">koch</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=chebsnake2">chebsnake2</a></b><b style="font-size:100%"> </b><b style="font-size:104%"><a href="/chebfun/examples/tags.php?query=autonomous">autonomous</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=Rosenbrock">Rosenbrock</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=ellipse">ellipse</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=least-squares">least-squares</a></b><b style="font-size:100%"> </b><b style="font-size:136%"><a href="/chebfun/examples/tags.php?query=heart">heart</a></b><b style="font-size:100%"> </b><b style="font-size:190%"><a href="/chebfun/examples/tags.php?query=gift">gift</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=Approximation">Approximation</a></b><b style="font-size:100%"> </b><b style="font-size:104%"><a href="/chebfun/examples/tags.php?query=chebfun2">chebfun2</a></b><b style="font-size:100%"> </b><b style="font-size:104%"><a href="/chebfun/examples/tags.php?query=advectiondiffusion">advectiondiffusion</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=encryption">encryption</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=Uniform distribution">Uniform distribution</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=Lebesgue">Lebesgue</a></b><b style="font-size:100%"> </b><b style="font-size:104%"><a href="/chebfun/examples/tags.php?query=blowup">blowup</a></b><b style="font-size:100%"> </b><b style="font-size:136%"><a href="/chebfun/examples/tags.php?query=poles">poles</a></b><b style="font-size:100%"> </b><b style="font-size:238%"><a href="/chebfun/examples/tags.php?query=quadrature">quadrature</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=logistic">logistic</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=asymptotics">asymptotics</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=singularities">singularities</a></b><b style="font-size:100%"> </b><b style="font-size:104%"><a href="/chebfun/examples/tags.php?query=orthogonal">orthogonal</a></b><b style="font-size:100%"> </b><b style="font-size:136%"><a href="/chebfun/examples/tags.php?query=parameter">parameter</a></b><b style="font-size:100%"> </b><b style="font-size:104%"><a href="/chebfun/examples/tags.php?query=calculus">calculus</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=game">game</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=gradient theorem">gradient theorem</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=infinite series">infinite series</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=extrapolation">extrapolation</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=Halphen">Halphen</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=9.28903">9.28903</a></b><b style="font-size:100%"> </b><b style="font-size:104%"><a href="/chebfun/examples/tags.php?query=boundarylayer">boundarylayer</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=greedy">greedy</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=mathjax">mathjax</a></b><b style="font-size:100%"> </b><b style="font-size:104%"><a href="/chebfun/examples/tags.php?query=Weierstrass">Weierstrass</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=amplitudemodulation">amplitudemodulation</a></b><b style="font-size:100%"> </b><b style="font-size:104%"><a href="/chebfun/examples/tags.php?query=rational function">rational function</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=message">message</a></b><b style="font-size:100%"> </b></span>
<h2 class="bigger" style="margin-left:-15pt; margin-top:10pt; margin-bottom:-10pt;">Function cloud</h2><br/>

<span style="width:400px"><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=SUM">SUM</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=INTERP1">INTERP1</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=3D">3D</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=LAGPOLY">LAGPOLY</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=ODE">ODE</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=LEBESGUE">LEBESGUE</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=PADEAPPROX">PADEAPPROX</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=PDE15S">PDE15S</a></b><b style="font-size:100%"> </b><b style="font-size:104%"><a href="/chebfun/examples/tags.php?query=REMEZ">REMEZ</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=SVD">SVD</a></b><b style="font-size:100%"> </b><b style="font-size:164%"><a href="/chebfun/examples/tags.php?query=CHEBPOLY">CHEBPOLY</a></b><b style="font-size:100%"> </b><b style="font-size:136%"><a href="/chebfun/examples/tags.php?query=RATINTERP">RATINTERP</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=AD">AD</a></b><b style="font-size:100%"> </b><b style="font-size:136%"><a href="/chebfun/examples/tags.php?query=ABS">ABS</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=CUMSUM">CUMSUM</a></b><b style="font-size:100%"> </b><b style="font-size:190%"><a href="/chebfun/examples/tags.php?query=DIRAC">DIRAC</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=SCRIBBLE">SCRIBBLE</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=CHEBELLIPSEPLOT">CHEBELLIPSEPLOT</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=CHEBPOLYPLOT">CHEBPOLYPLOT</a></b><b style="font-size:100%"> </b><b style="font-size:136%"><a href="/chebfun/examples/tags.php?query=CONV">CONV</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=EQUI">EQUI</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=HERMPOLY">HERMPOLY</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=LEGPOLY">LEGPOLY</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=PLOT3">PLOT3</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=ARCLENGTH">ARCLENGTH</a></b><b style="font-size:100%"> </b><b style="font-size:66%"><a href="/chebfun/examples/tags.php?query=AM">AM</a></b><b style="font-size:100%"> </b></span>
</div> <!--END CONTENT-->
<div class="hr1"> <hr /> </div>

	<div id="myfooter"> <!-- FOOTER -->
		<p>Please <a href="/chebfun/contact/">contact us</a> with any questions and comments.<br/>
		Copyright &copy; 2013, The University of Oxford &amp; The Chebfun Team.</p>
<!--		<p style="color: #AAA;"><a href="http://validator.w3.org/check?uri=referer" style="color: #888;">XHTML 1.0 Transitional</a> and <a href="http://jigsaw.w3.org/css-validator/check/referer" style="color: #888;">CSS 2.1</a> compliant.</p>
-->

	</div> <!-- END FOOTER -->
</div> <!--END WRAPPER-->
<div id="footimg"></div>
</body>
</html>