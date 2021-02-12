////////////////////////////////////////////////////////////////////
//   SCRIPT TO SELECT, DISPLAY AND EXPORT SENTINEL-2 SWIR DATA    //
//    (c) B. Smets, Royal Museum for Central Africa, Belgium      //
//                                                                //
// CITATION:                                                      //
// Smets, B. (2020) - Select, Display and Export Sentinel-2 SWIR  //
// image over a volcano, in Google Earth Engine.                  //
// https://github.com/besmets/GEE_Scripts.                        //
// Accessed online on DD/MM/YYYY.                                 //
//                                                                //
// Last modified: 30th March 2020                                 //
////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////
//                     VOLCANO = NYAMULAGIRA                      //
////////////////////////////////////////////////////////////////////

//Define the dates of Interest
//------------------------------------------------------------------
var start = ee.Date("2020-01-01");   ////////// TO MODIFY //////////
var finish = ee.Date("2020-12-31");  ////////// TO MODIFY //////////

//Define the image to display
//------------------------------------------------------------------
var image_select_last = 'manual';      /// write 'auto' or 'manual' //
var image_number = 0;  //manual selection /////// TO MODIFY ///////

//Exporting option
//------------------------------------------------------------------
var roi_export = 'no';              ////////// TO MODIFY //////////

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

// Define the Region of Interest (ROI)
//------------------------------------------------------------------
var roi = ee.Geometry.Rectangle([29.165, -1.37, 29.245, -1.45]);

// Define the location of the volcano and the zoom
//------------------------------------------------------------------
var locx = 29.2;
var locy = -1.41;
var zoom = 14;

// Filter the Sentinel-2 collection
//------------------------------------------------------------------

var collection = "COPERNICUS/S2";
var s2_coll = ee.ImageCollection(collection)
.filterDate(start, finish)
.filterBounds(roi);

// List the collection and display it in the console
var imagelist = s2_coll.toList(s2_coll.size());
print('List of images: ', imagelist);

// Automatically select the last image of the list
var list_length_raw = ee.String(imagelist.length());   // string of the image length
var list_length_str = list_length_raw.getInfo();       // extract the string from GEE
var list_length_int = parseInt(list_length_str, 10);   // string to integer (10 = decimal radix)
var image_number_auto = list_length_int-1;

// List the image IDs and display it in the console
var idlist = imagelist.map(function(item) {
  return ee.Image(item).id();
});
print('List of IDs: ', idlist);

// Image Display
//------------------------------------------------------------------

//Provide image ID
if (image_select_last == 'auto') {
  var imageID1 = imagelist.get(image_number_auto);
  print('Selected image (auto) =', imageID1);
}
else {
  var imageID1 = imagelist.get(image_number);
  print('Selected image =', imageID1);
}

// Load an image.
var image1 = ee.Image(imageID1);

//Extract the date of the selected image in YYYY-MM-DD format
var image_date2 = ee.Date(image1.get('system:time_start')).format('YYYY-MM-dd');

// Define the visualization parameters.
var visParams = {
  bands: ['B12', 'B11', 'B8A'],
  min: -100,
  max: 2500,
};

// Center the map to the volcano and display the cropped image
//------------------------------------------------------------------
Map.setCenter(locx, locy, zoom);
Map.addLayer(image1.clip(roi), visParams);

// IMAGE EXPORT OPTION
//------------------------------------------------------------------

if (roi_export == 'yes') {
  var description = 'Export_SWIR_image_' + image_date2.getInfo();
  var export_name = 'S2_Nyam_B12118A_' + image_date2.getInfo();
  var image_to_export = image1.select(['B12', 'B11', 'B8A'])
  Export.image.toDrive({
    image: image_to_export,
    description: description,
    folder: '_fromGEE/_S2_snapshot_NYAM',
    fileNamePrefix: export_name,
    scale: 20,
    region: roi,
    fileFormat: 'GeoTIFF'
  });
}