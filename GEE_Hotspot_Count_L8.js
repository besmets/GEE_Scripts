/////////////////////////////////////////////////////////////////////
//       SCRIPT TO COUNT HOTSPOT PIXELS IN LANDSAT-8 IMAGES        //
//    (c) B. Smets, Royal Museum for Central Africa, Belgium       //
//                                                                 //
// CITATION:                                                       //
// Smets, B. (2021) - Script to count hotspot pixels in            //
// Landsat-8 images, in Google Earth Engine.                       //
// https://github.com/besmets/GEE_Scripts.                         //
// Accessed online on DD/MM/YYYY.                                  //
//                                                                 //
// Last modified: 20th January 2021                                //
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
//                       VOLCANO = KILAUEA                         //
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
///////////////////////////   Variables    //////////////////////////
/////////////////////////////////////////////////////////////////////

//Define the volcano and the region of interest (roi)
var vulc='KILAUEA (SUMMIT)'; //COMMENT OR DECOMMENT
var roi = ee.Geometry.Rectangle([-155.31, 19.38, -155.25, 19.44]); // Kīlauea summit

//Define the dates of Interest
var start = ee.Date("2020-06-01");
var finish = ee.Date("2021-01-19");

/////////////////////////////////////////////////////////////////////
///////////////////////// SELECT COLLECTION /////////////////////////
/////////////////////////////////////////////////////////////////////

print('## ANALYZING HOTSPOTS ##');
print('##   FOR ' + vulc + '  ##');

// Define the location of the volcano
var loc = ee.Geometry.Point([-155.28, 19.41]);
 
/// THE LANDSAT-8 COLLECTION ///
var collection = "LANDSAT/LC08/C01/T1_TOA";
var l8_coll = ee.ImageCollection(collection)
.filterDate(start, finish)
.filterBounds(loc);
 
// List the collection and display it
var imagelist = l8_coll.toList(l8_coll.size());
print('List of images: ', imagelist);
 
// List the image IDs and display it
var idlist = imagelist.map(function(item) {
  return ee.Image(item).id();
});
print('List of IDs: ', idlist);

//List of dates and display it
var datelist= imagelist.map(function(item){
  var id=ee.Image(item).id();
  var year=id.slice(12,16);
  var month=id.slice(16,18);
  var day=id.slice(18,20);
  return ee.Feature(ee.Geometry.Point(0,0),{value: ee.String(year).cat(month).cat(day)});
});
print('List of dates: ',datelist);
var FCdates = ee.FeatureCollection(ee.List(datelist));

//Initialize list
var LIST = ee.List([]);
 
/////////////////////////////////////////////////////////////////////
/////////////////////////// THE FUNCTION // /////////////////////////
/////////////////////////////////////////////////////////////////////

function HTDfunction(image,LIST){
  
  ////////////////////// CALCULATE ALPHA AND BETA //////////////////////
  
  //Define the conditions for alpha
  var alpha1 = image.expression('B7/B6', {'B7': image.select('B7'),'B6': image.select('B6')});
  var alpha2 = image.expression('B7/B5', {'B7': image.select('B7'),'B5': image.select('B5')});
  
  //Define the conditions for beta
  var beta1 = image.expression('B6/B5', {'B6': image.select('B6'),'B5': image.select('B5')});
  
  //Select the bands needed for the calculation of alpha and beta
  var B7 = image.select('B7');
  var B6 = image.select('B6');
  
  //Calculate alpha based on the conditions
  var alpha = image.where(alpha1.gte(1.4), 1);
  var alpha = alpha.where(alpha.neq(1),0);
  
  var alphaa = image.where(alpha2.gte(1.4),1);
  var alphaa = alphaa.where(alphaa.neq(1),0);
  
  var alphaaa = image.where(B7.gte(0.15),1);
  var alphaaa = alphaaa.where(alphaaa.neq(1),0);
  
  var Alphaa = image.expression(
      'alpha+alphaa+alphaaa', {
        'alpha': alpha,
        'alphaa': alphaa,
        'alphaaa': alphaaa
  });
  var Alpha = Alphaa.where(Alphaa.eq(3),1);
  var Alpha = Alpha.where(Alphaa.neq(3),0);
  
  //Calculate beta based on the conditions
  var s = image.where(B6.lt(0),1);
  var ss = s.where(B7.lt(0),1);
  var S = ss.where(ss.neq(1),0);
  
  var beta = image.where(beta1.gte(2),1);
  var beta = beta.where(beta.neq(1),0); 
  
  var betaa = image.where(B6.gte(0.5),1);
  var betaa = betaa.where(betaa.neq(1),0);
  
  var sumbetaa =  image.expression(
      'beta+betaa', {
        'beta': beta,
        'betaa': betaa,
  });
  var betaaa = sumbetaa.where(sumbetaa.eq(2),1);
  var betaaa = betaaa.where(betaaa.neq(1),0);
  
  var betaaaa = betaaa.where (betaaa.eq(1),1);
  var betaaaa =  betaaaa.where(S.eq(1),1);
  
  var Beta = betaaaa.where(betaaaa.neq(1),0);

  ////////////////////////////// CLUSTER //////////////////////////////
  
  // Create a list of weights for a 3x3 kernel.
  var list = [1, 1, 1];
  // The center of the kernel is zero.
  var centerList = [1, 0, 1];
  // Assemble a list of lists: the 3x3 kernel weights as a 2-D matrix.
  var lists = [list, centerList, list];
  // Create the kernel from the weights.
  var kernel = ee.Kernel.fixed(3, 3, lists, -4, -4, false);
  
  var cluster= Alpha.convolve(kernel);
  var Cluster = cluster.where(Beta.eq(0),-100000);
  
  var hotspot = image.where(Alpha.eq(0),0);
  var hotspot1 = hotspot.where(Beta.eq(0),0);
  var hotspot2 = hotspot1.where(Alpha.eq(1),1);
  var Hotspot = hotspot2.where(Cluster.gte(1),1);
  
  ////////////////////////////// MASK //////////////////////////////

  var clipped= Hotspot.clip(roi);
  var Mask= clipped.selfMask();
  
  /////////////////////////////// COUNT ///////////////////////////////
  //improved count
  var count = Mask.reduceRegion({
    reducer: ee.Reducer.count(),
    geometry: Mask.geometry(),
    scale: 30
  });
  return ee.List(LIST).add(ee.Feature(ee.Geometry.Point(0,0),{value: count.get('B1')}));
}

/////////////////////////////////////////////////////////////////////
//////////////// CALCULATE PIXEL COUNT WITH FUNCTION ////////////////
/////////////////////////////////////////////////////////////////////

var maskcollection=l8_coll.iterate(HTDfunction,LIST);
maskcollection = ee.FeatureCollection(ee.List(maskcollection));

/////////////////////////////////////////////////////////////////////
/////////////////// JOIN TABLES AND EXPORT TO CSV ///////////////////
/////////////////////////////////////////////////////////////////////

//join tables
// Use an equals filter to define how the collections match.
var filter = ee.Filter.equals({
  leftField: 'system:index',
  rightField: 'system:index'
});

// Create the join.
var innerJoin = ee.Join.inner();

// Apply the join.
var Joined = innerJoin.apply( FCdates, maskcollection, filter);
print (Joined);

//export to a csv
var str_start = start.format('YYYYMMdd').getInfo();
var str_finish = finish.format('YYYYMMdd').getInfo();
var filename = 'L8_hotspot_count_Kilauea_summit_'+str_start+'_'+str_finish;
Export.table.toDrive({
  collection: Joined,
  description: filename,
  folder: 'fromGEE',
  fileNamePrefix: filename,
  fileFormat: 'CSV',
  selectors: ['primary.value, secondary.value']
});
