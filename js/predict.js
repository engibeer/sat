// Globals

var km2mi = 0.621371;
var deg2rad = 1.745329251994330e-2;
var pi = 3.14159265358979323846;
var pio2 = 1.57079632679489656;
var x3pio2 = 4.71238898038468967;
var twopi = 6.28318530717958623;
var e6a = 1.0e-6;
var twothrd = 6.6666666666666666e-1;
var xj2 = 1.0826158e-3;
var xj3 = -2.53881e-6;
var xj4 = -1.65597e-6;
var xke = 7.43669161e-2;
var xkmper = 6.378137e3;
var xmnpda = 1.44e3;
var ae = 1.0;
var ck2 = 5.413079e-4;
var ck4 = 6.209887e-7;
var f = 3.35281066474748e-3;
var ge = 3.986008e5;
var s = 1.012229;
var qoms2t = 1.880279e-09;
var secday = 8.6400e4;
var omega_e = 1.00273790934;
var omega_eR = 6.3003879;
var zns = 1.19459e-5;
var c1ss = 2.9864797e-6;
var zes = 1.675e-2;
var znl = 1.5835218e-4;
var c1l = 4.7968065e-7;
var zel = 5.490e-2;
var zcosis = 9.1744867e-1;
var zsinis = 3.9785416e-1;
var zsings = -9.8088458e-1;
var zcosgs = 1.945905e-1;
var zcoshs = 1;
var zsinhs = 0;
var q22 = 1.7891679e-6;
var q31 = 2.1460748e-6;
var q33 = 2.2123015e-7;
var g22 = 5.7686396;
var g32 = 9.5240898e-1;
var g44 = 1.8014998;
var g52 = 1.0508330;
var g54 = 4.4108898;
var root22 = 1.7891679e-6;
var root32 = 3.7393792e-7;
var root44 = 7.3636953e-9;
var root52 = 1.1428639e-7;
var root54 = 2.1765803e-9;
var thdt = 4.3752691e-3;
var rho = 1.5696615e-1;
var mfactor = 7.292115e-5;
var sr = 6.96000e5;
var AU = 1.49597870691e8;

function TwoLineElement(name, line1, line2){
	// see http://en.wikipedia.org/wiki/Two-line_elements
	this.name = name;
	this.line1 = line1;
	this.line2 = line2;
	this.sat_num = line1.substring(2,7);
	this.classification = line1.substring(7,8);
	this.launch_year = line1.substring(9,11);
	this.launch_number = line1.substring(11,14);
	this.launch_piece = line1.substring(14,17);
	this.epoch_year = line1.substring(18,20);
	this.epoch = line1.substring(20,32);
	this.dt_mm_2 = line1.substring(33,43);
	this.d2t_mm_6 = 1.0e-5 * line1.substring(45,50) / Math.pow(10.0, line1.substring(50,52));
	this.bstar_drag = 1.0e-5 * line1.substring(53,59) / Math.pow(10.0, line1.substring(60, 61));
	this.inclination = line2.substring(8,16);
	this.right_ascension = line2.substring(17,25);
	this.eccentricity = 1.0e-07 * line2.substring(26,33);
	this.arg_perigee = line2.substring(34,42);
	this.mean_anomaly = line2.substring(43,51);
	this.mean_motion = line2.substring(52,63);
	this.orbitnum = line2.substring(63,68);

	// degrees to radians
	this.right_ascension *= deg2rad;
	this.arg_perigee *= deg2rad;
	this.mean_anomaly *= deg2rad;
	this.inclination *= deg2rad;
	var temp = (2*pi/xmnpda) / xmnpda;
	this.mean_motion = this.mean_motion * (2*pi/xmnpda);
	this.dt_mm_2 *= temp;
	this.d2t_mm_6 *= (2*pi/xmnpda);
	this.period = 1440.0 / this.mean_motion;
	var sma = 331.25 * Math.exp(Math.log(1440.0 / this.mean_motion) * (2.0/3.0));
	var c1 = Math.cos(this.inclination * deg2rad);
	var e2 = 1.0 - (this.eccentricity * this.eccentricity);
	this.nodal_period = (this.period * 360.0) / 
						(360.0 + (4.97 * Math.pow(xkmper / sma, 3.5) * (5.0 * c1 * c1 - 1.0) / (e2 * e2) / this.mean_motion));

	// period > 225 is a deep space orbit
	var dd1 = xke / this.mean_motion;
	var a1 = Math.pow(dd1, (2.0 / 3.0));
	var r1 = Math.cos(this.inclination);
	dd1 = 1.0 - this.eccentricity * this.eccentricity;
	temp = ck2 * 1.5 * (r1 * r1 * 3.0 - 1.0) / Math.pow(dd1, 1.5);
	var del1 = temp / (a1 * a1);
	var ao = a1 * (1.0 - del1 * 2.0/3.0 * .5) + del1 * (del1 * 1.654320987654321 + 1.0);
	var delo = temp / (ao*ao);
	var xnodp = this.mean_motion / (delo + 1.0);

	// deep space or near earth ephemeris
	if(2*pi / xnodp / xmnpda >= .15625)
		this.deep = true;
	else
		this.deep = false;
}
// geodetic object
function Geodetic(){
	this.lat = 0;
	this.long = 0;
	this.alt = 0;
	this.theta = 0;
}

// constructor
function Predict(){
	
}

Predict.prototype.AosHappens = function(tle, obs_geodetic){
	if(tle.mean_motion == 0.0)
		return false;
	else
		return true;
}

Predict.prototype.Geostatioary = function(tle){
	// TODO
}

Predict.prototype.Decayed = function(tle){
	// TODO
}

// This function finds and returns the time of AOS (aostime).
Predict.prototype.FindAos = function(tle){
	// TODO
}

// returns next AOS time for given TLE and Geodetic position
Predict.prototype.NextAos = function(tle){
	// TODO
}

Predict.prototype.SGP4 = function(tsince, tle){
	// tsince: time since epoch in minutes
	// return pos and vel vectors as ECI satellite position and velocity

	var aodp;
	var aycof;
	var c1;
	var c4;
	var c5;
	var cosio;
	var d2;
	var d3;
	var d4;
	var delmo;
	var omgcof;
	var eta;
	var omgdot;
	var sinio;
	var xnodp;
	var sinmo;
	var t2cof;
	var t3cof;
	var t4cof;
	var t5cof;
	var x1mth2;
	var x3thm1;
	var x7thm1;
	var xmcof;
	var xmdot;
	var xnodcf;
	var xnodot;
	var xlcof;

	// Initialize
	// Recover original mean motion (xnodp) and
	// semimajor axis (aodp) from input elements.
	a1 = Math.pow(xke / tle.mean_motion, (2.0/3.0));
	cosio = Math.cos(tle.inclination);
	theta2 = cosio * cosio;
	x3thm1 = 3 * theta2 - 1.0;
	eosq = tle.eccentricity * tle.eccentricity;
	betao2 = 1.0 - eosq;
	betao = Math.sqrt(betao2);
	del1 = 1.5 * ck2 * x3thm1 / (a1 * a1 * betao * betao2);
	ao = a1 * (1.0 - del1 * (0.5 * (2.0/3.0) + del1 * (1.0 + 134.0 / 81.0 * del1)));
	delo = 1.5 * ck2 * x3thm1 / (ao * ao * betao * betao2);
	xnodp = tle.mean_motion / (1.0 + delo);
	aodp = ao / (1.0 - delo);


	// For perigee less than 220 kilometers, the "simple"
	// flag is set and the equations are truncated to linear
	// variation in sqrt a and quadratic variation in mean
	// anomaly.  Also, the c3 term, the delta omega term, and
	// the delta m term are dropped.

	if((aodp * (1 - tle.eccentricity) / ae) < (220 / xkmper + ae))
		this.SIMPLE_FLAG = 1;
	else
		this.SIMPLE_FLAG = 0;

	// For perigees below 156 km, the
	// values of s and qoms2t are altered.
	var s4 = s;
	var qoms24 = qoms2t
	var perigee = (aodp * (1 - tle.eccentricity) - ae) * xkmper;
	if(perigee < 156.0) {
		if(perigee <= 98.0)
			s4 = 20;
		else
			s4 = perigee - 78.0;
		qoms24 = Math.pow((120 - s4) * ae / xkmper, 4);
		s4 = s4 / xkmper + ae;
	}

	// shit is about to get real
	pinvsq = 1 / (aodp * aodp * betao2 * betao2);
	tsi = 1 / (aodp - s4);
	eta = aodp * tle.eccentricity * tsi;
	etasq = eta * eta;
	eeta = tle.eccentricity * eta;
	psisq = Math.abs(1 - etasq);
	coef = qoms24 * Math.pow(tsi, 4);
	coef1 = coef / Math.pow(psisq, 3.5);
	c2 = coef1 * xnodp * (aodp * (1 + 1.5 * etasq + eeta * (4 + etasq)) + 0.75 * ck2 * tsi / psisq * x3thm1 * (8 + 3 * etasq * (8 + etasq)));
	c1 = tle.bstar_drag * c2;
	sinio = Math.sin(tle.inclination);
	a3ovk2 = -xj3 / ck2 * Math.pow(ae, 3);
	c3 = coef * tsi * a3ovk2 * xnodp * ae * sinio / tle.eccentricity;
	x1mth2 = 1 - theta2;

	c4 = 2 * xnodp * coef1 * aodp * betao2 * (eta * (2 + 0.5 * etasq) + tle.eccentricity * (0.5 + 2 * etasq) - 2 * ck2 * tsi / (aodp * psisq) * (-3 * x3thm1 * (1 - 2 * eeta + etasq * (1.5 - 0.5 * eeta)) + 0.75 * x1mth2 * (2 * etasq - eeta * (1 + etasq)) * Math.cos(2 * tle.arg_perigee)));
	c5 = 2 * coef1 * aodp * betao2 * (1 + 2.75 * (etasq + eeta) + eeta * etasq);
	theta4 = theta2 * theta2;
	temp1 = 3 * ck2 * pinvsq * xnodp;
	temp2 = temp1 * ck2 * pinvsq;
	temp3 = 1.25 * ck4 * pinvsq * pinvsq * xnodp;
	xmdot = xnodp + 0.5 * temp1 * betao * x3thm1 + 0.0625 * temp2 * betao * (13 - 78 * theta2 + 137 * theta4);
	x1m5th = 1 - 5 * theta2;
	omgdot = -0.5 * temp1 * x1m5th + 0.0625 * temp2 * (7 - 114 * theta2 + 395 * theta4) + temp3 * (3 - 36 * theta2 + 49 * theta4);
	xhdot1 = -temp1 * cosio;
	xnodot = xhdot1 + (0.5 * temp2 * (4 - 19 * theta2) + 2 * temp3 * (3 - 7 * theta2)) * cosio;
	omgcof = tle.bstar_drag * c3 * Math.cos(tle.arg_perigee);
	xmcof = -(2.0/3.0) * coef * tle.bstar_drag * ae / eeta;
	xnodcf = 3.5 * betao2 * xhdot1 * c1;
	t2cof = 1.5 * c1;
	xlcof = 0.125 * a3ovk2 * sinio * (3 + 5 * cosio) / (1 + cosio);
	aycof = 0.25 * a3ovk2 * sinio;
	delmo = Math.pow(1 + eta * Math.cos(tle.mean_anomaly), 3);
	sinmo = Math.sin(tle.mean_anomaly);
	x7thm1 = 7 * theta2 - 1;

	if(!this.SIMPLE_FLAG) {
		c1sq = c1 * c1;
		d2 = 4 * aodp * tsi * c1sq;
		temp = d2 * tsi * c1 / 3;
		d3 = (17 * aodp + s4) * temp;
		d4 = 0.5 * temp * aodp * tsi * (221 * aodp + 31 * s4) * c1;
		t3cof = d2 + 2 * c1sq;
		t4cof = 0.25 * (3 * d3 + c1 * (12 * d2 + 10 * c1sq));
		t5cof = 0.2 * (3 * d4 + 12 * c1 * d3 + 6 * d2 * d2 + 15 * c1sq * (2 * d2 + c1sq));
	}

	// Update for secular gravity and atmospheric drag.

	xmdf = tle.mean_anomaly + xmdot * tsince;
	omgadf = tle.arg_perigee + omgdot * tsince;
	xnoddf = tle.right_ascension + xnodot * tsince;
	omega = omgadf;
	xmp = xmdf;
	tsq = tsince * tsince;
	xnode = xnoddf + xnodcf * tsq;
	tempa = 1 - c1 * tsince;
	tempe = tle.bstar_drag * c4 * tsince;
	templ = t2cof * tsq;

	if(!this.SIMPLE_FLAG) {
		delomg = omgcof * tsince;
		delm = xmcof * (Math.pow(1 + eta * Math.cos(xmdf), 3) - delmo);
		temp = delomg + delm;
		xmp = xmdf + temp;
		omega = omgadf - temp;
		tcube = tsq * tsince;
		tfour = tsince * tcube;
		tempa = tempa - d2 * tsq - d3 * tcube - d4 * tfour;
		tempe = tempe + tle.bstar_drag * c5 * (Math.sin(xmp) - sinmo);
		templ = templ + t3cof * tcube + tfour * (t4cof + tsince * t5cof);
	}

	a = aodp * Math.pow(tempa, 2);
	e = tle.eccentricity - tempe;
	xl = xmp + omega + xnode + xnodp * templ;
	beta = Math.sqrt(1 - e * e);
	xn = xke / Math.pow(a, 1.5);

	// Long period periodics

	axn = e * Math.cos(omega);
	temp = 1 / (a * beta * beta);
	xll = temp * xlcof * axn;
	aynl = temp * aycof;
	xlt = xl + xll;
	ayn = e * Math.sin(omega) + aynl;

	// Solve Kepler's Equation
	// http://mathworld.wolfram.com/KeplersEquation.html

	capu = (xlt - xnode) % (2*Math.PI);
	temp2 = capu;
	i = 0;

	do {
		sinepw = Math.sin(temp2);
		cosepw = Math.cos(temp2);
		temp3 = axn * sinepw;
		temp4 = ayn * cosepw;
		temp5 = axn * cosepw;
		temp6 = ayn * sinepw;
		epw = (capu - temp4 + temp3 - temp2) / (1 - temp5 - temp6) + temp2;

		if(Math.abs(epw - temp2) <= e6a)
			break;

		temp2 = epw;
	} while (i++<10);

	// Short period preliminary quantities

	ecose = temp5 + temp6;
	esine = temp3 - temp4;
	elsq = axn * axn + ayn * ayn;
	temp = 1 - elsq;
	pl = a * temp;
	r = a * (1 - ecose);
	temp1 = 1 / r;
	rdot = xke * Math.sqrt(a) * esine * temp1;
	rfdot = xke * Math.sqrt(pl) * temp1;
	temp2 = a * temp1;
	betal = Math.sqrt(temp);
	temp3 = 1 / (1 + betal);
	cosu = temp2 * (cosepw - axn + ayn * esine * temp3);
	sinu = temp2 * (sinepw - ayn - axn * esine * temp3);
	u = AcTan(sinu, cosu);
	sin2u = 2 * sinu * cosu;
	cos2u = 2 * cosu * cosu - 1;
	temp = 1 / pl;
	temp1 = ck2 * temp;
	temp2 = temp1 * temp;

	// Update for short periodics

	rk = r * (1 - 1.5 * temp2 * betal * x3thm1) + 0.5 * temp1 * x1mth2 * cos2u;
	uk = u - 0.25 * temp2 * x7thm1 * sin2u;
	xnodek = xnode + 1.5 * temp2 * cosio * sin2u;
	xinck = tle.inclination + 1.5 * temp2 * cosio * sinio * cos2u;
	rdotk = rdot - xn * temp1 * x1mth2 * sin2u;
	rfdotk = rfdot + xn * temp1 * (x1mth2 * cos2u + 1.5 * x3thm1);

	// Orientation vectors

	sinuk = Math.sin(uk);
	cosuk = Math.cos(uk);
	sinik = Math.sin(xinck);
	cosik = Math.cos(xinck);
	sinnok = Math.sin(xnodek);
	cosnok = Math.cos(xnodek);
	xmx = -sinnok * cosik;
	xmy = cosnok * cosik;
	ux = xmx * sinuk + cosnok * cosuk;
	uy = xmy * sinuk + sinnok * cosuk;
	uz = sinik * sinuk;
	vx = xmx * cosuk - cosnok * sinuk;
	vy = xmy * cosuk - sinnok * sinuk;
	vz = sinik * cosuk;

	// Position and velocity

	var x_pos = rk * ux;
	var y_pos = rk * uy;
	var z_pos = rk * uz;

	var x_vel = rdotk * ux + rfdotk * vx;
	var y_vel = rdotk * uy + rfdotk * vy;
	var z_vel = rdotk * uz + rfdotk * vz;

	var result = {
		pos: { x: x_pos, y: y_pos, z: z_pos },
		vel: { x: x_vel, y: y_vel, z: z_vel }
	};

	return result;
}

// Four-quadrant arctan function
function AcTan(sinx, cosx) {
	if(cosx == 0.0) {
		if(sinx > 0.0)
			return (pio2);
		else
			return (x3pio2);
	} else {
		if(cosx > 0.0) {
			if(sinx > 0.0)
				return (Math.atan(sinx / cosx));
			else
				return (twopi + Math.atan(sinx / cosx));
		} else
			return (pi + Math.atan(sinx / cosx));
	}
}

// returns minutes since epoch of given tle
function since_epoch(tle){
	var year = tle.epoch_year;
	var days = tle.epoch;
	var today = new Date();
	if(year < 56)
		year = parseInt(year) + 2000;
	else
		year = parseInt(year) + 1900;
	var year_diff = today.getFullYear() - year;
	var minutes = (year_diff * 365 * 24 * 60) + today.getMinuteOfYear() - parseFloat(days);
	return minutes;
}

function Convert_Sat_State(pos, vel) {
	// Converts the satellite's position and velocity
	// vectors from normalized values to km and km/sec
	var s_pos = Scale_Vector(xkmper, pos);
	var s_vel = Scale_Vector(xkmper * xmnpda / secday, vel);
	return [s_pos, s_vel];
}

function Scale_Vector(k, v) {
	// Multiplies the vector v1 by the scalar k
	v.x *= k;
	v.y *= k;
	v.z *= k;
	return Magnitude(v);
}

function Magnitude(v) {
	// Calculates scalar magnitude of a vector argument
	return Math.sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
}

Date.prototype.getMinuteOfYear = function() {
	var onejan = new Date(this.getFullYear(),0,1);
	return (this - onejan) / 60000;
}

function calc_ISS(){
	var name = "ISS (ZARYA)";
	var line1 = "1 25544U 98067A   11331.86429521  .00069397  00000-0  86074-3 0  5163"
	var line2 = "2 25544  51.6420  57.1271 0030727  98.4286   5.8951 15.58813992746465";
	tle = new TwoLineElement(name, line1, line2);
	predict = new Predict(tle);
	tsince = since_epoch(tle);
	result = predict.SGP4(tsince, tle);
	// note js passes objects by reference,
	// this clobbers original result
	Convert_Sat_State(result.pos, result.vel);
	return result;
}