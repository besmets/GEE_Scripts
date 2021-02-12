/////////////////////////////////////////////////////////////////////
//      SCRIPT TO COUNT HOTSPOT PIXELS IN SENTINEL-2 IMAGES        //
//    (c) B. Smets, Royal Museum for Central Africa, Belgium       //
//                                                                 //
// CITATION:                                                       //
// Smets, B. (2021) - Script to count hotspot pixels in            //
// Sentinel-2 images, in Google Earth Engine.                      //
// https://github.com/besmets/GEE_Scripts.                         //
// Accessed online on DD/MM/YYYY.                                  //
//                                                                 //
// Last modified: 20th january 2021                                //
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
//                       VOLCANO = KILAUEA                         //
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
///////////////////////////   Variables    //////////////////////////
/////////////////////////////////////////////////////////////////////

//Define the volcano and the region of interest (roi)
var vulc='KILAUEA (SUMMIT)'; //COMMENT OR DECOMMENT
var roi = ee.Geometry.Rectangle([-155.31, 19.38, -155.25, 19.44]); // Kīlauea

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
 
/// THE SENTINEL-2 COLLECTION ///
var collection = "COPERNICUS/S2";
var s2_coll = ee.ImageCollection(collection)
.filterDate(start, finish)
.filterBounds(loc);
 
// List the collection and display it
var imagelist = s2_coll.toList(s2_coll.size());
print('List of images: ', imagelist);
 
// List the image IDs and display it
var idlist = imagelist.map(function(item) {
  return ee.Image(item).id();
});
print('List of IDs: ', idlist);

//List of dates and display it
var datelist= imagelist.map(function(item){
  var id=ee.Image(item).id();
  var year=id.slice(0,4);
  var month=id.slice(4,6);
  var day=id.slice(6,8);
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
  var alpha1 = image.expression('B12/B11', {'B12': image.select('B12'),'B11': image.select('B11')});
  var alpha2 = image.expression('B12/B8A', {'B12': image.select('B12'),'B8A': image.select('B8A')});
  
  //Define the conditions for beta
  var beta1 = image.expression('B11/B8A', {'B11': image.select('B11'),'B8A': image.select('B8A')});
  
  //Select the bands needed for the calculation of alpha and beta
  var B12 = image.select('B12');
  var B11 = image.select('B11');
  
  //Calculate alpha based on the conditions
  var alpha_A_raw = image.where(alpha1.gte(1.4), 1);
  var alpha_A = alpha_A_raw.where(alpha_A_raw.neq(1),0);
  
  var alpha_B_raw = image.where(alpha2.gte(1.4),1);
  var alpha_B = alpha_B_raw.where(alpha_B_raw.neq(1),0);
  
  var alpha_C_raw = image.where(B12.gte(0.15),1);
  var alpha_C = alpha_C_raw.where(alpha_C_raw.neq(1),0);
  
  var Alpha_total = image.expression(
      'alpha+alphaa+alphaaa', {
        'alpha': alpha_A,
        'alphaa': alpha_B,
        'alphaaa': alpha_C
  });
  var Alpha_raw = Alpha_total.where(Alpha_total.eq(3),1);
  var Alpha = Alpha_raw.where(Alpha_total.neq(3),0);
  
  //Calculate beta based on the conditions
  var s = image.where(B11.lt(0),1);
  var ss = s.where(B12.lt(0),1);
  var S = ss.where(ss.neq(1),0);
  
  var beta_A_raw = image.where(beta1.gte(2),1);
  var beta_A = beta_A_raw.where(beta_A_raw.neq(1),0); 
  
  var beta_B_raw = image.where(B11.gte(0.5),1);
  var beta_B = beta_B_raw.where(beta_B_raw.neq(1),0);
  
  var beta_AB =  image.expression( //previsously sumbetaa
      'beta+betaa', {
        'beta': beta_A,
        'betaa': beta_B,
  });
  var beta_C_raw = beta_AB.where(beta_AB.eq(2),1);
  var beta_C = beta_C_raw.where(beta_C_raw.neq(1),0);
  
  var beta_D_raw = beta_C.where (beta_C.eq(1),1);
  var beta_D =  beta_D_raw.where(S.eq(1),1);
  
  var Beta = beta_D.where(beta_D.neq(1),0);

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
    scale: 20
  });
  return ee.List(LIST).add(ee.Feature(ee.Geometry.Point(0,0),{value: count.get('B1')}));
}

/////////////////////////////////////////////////////////////////////
//////////////// CALCULATE PIXEL COUNT WITH FUNCTION ////////////////
/////////////////////////////////////////////////////////////////////

var maskcollection=s2_coll.iterate(HTDfunction,LIST);
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
var filename = 'S2_hotspot_count_Kilauea_summit_'+str_start+'_'+str_finish;
Export.table.toDrive({
  collection: Joined,
  description: filename,
  folder: 'fromGEE',
  fileNamePrefix: filename,
  fileFormat: 'CSV',
  selectors: ['primary.value, secondary.value']
});
