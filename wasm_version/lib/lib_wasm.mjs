
var wasm_ = (() => {
  var _scriptDir = import.meta.url;
  
  return (
function(wasm_) {
  wasm_ = wasm_ || {};

var Module=typeof wasm_!="undefined"?wasm_:{};var readyPromiseResolve,readyPromiseReject;Module["ready"]=new Promise(function(resolve,reject){readyPromiseResolve=resolve;readyPromiseReject=reject});var moduleOverrides=Object.assign({},Module);var arguments_=[];var thisProgram="./this.program";var quit_=(status,toThrow)=>{throw toThrow};var ENVIRONMENT_IS_WEB=typeof window=="object";var ENVIRONMENT_IS_WORKER=typeof importScripts=="function";var ENVIRONMENT_IS_NODE=typeof process=="object"&&typeof process.versions=="object"&&typeof process.versions.node=="string";var scriptDirectory="";function locateFile(path){if(Module["locateFile"]){return Module["locateFile"](path,scriptDirectory)}return scriptDirectory+path}var read_,readAsync,readBinary,setWindowTitle;function logExceptionOnExit(e){if(e instanceof ExitStatus)return;let toLog=e;err("exiting due to exception: "+toLog)}var fs;var nodePath;var requireNodeFS;if(ENVIRONMENT_IS_NODE){if(ENVIRONMENT_IS_WORKER){scriptDirectory=require("path").dirname(scriptDirectory)+"/"}else{scriptDirectory=__dirname+"/"}requireNodeFS=()=>{if(!nodePath){fs=require("fs");nodePath=require("path")}};read_=function shell_read(filename,binary){requireNodeFS();filename=nodePath["normalize"](filename);return fs.readFileSync(filename,binary?undefined:"utf8")};readBinary=filename=>{var ret=read_(filename,true);if(!ret.buffer){ret=new Uint8Array(ret)}return ret};readAsync=(filename,onload,onerror)=>{requireNodeFS();filename=nodePath["normalize"](filename);fs.readFile(filename,function(err,data){if(err)onerror(err);else onload(data.buffer)})};if(process["argv"].length>1){thisProgram=process["argv"][1].replace(/\\/g,"/")}arguments_=process["argv"].slice(2);process["on"]("uncaughtException",function(ex){if(!(ex instanceof ExitStatus)){throw ex}});process["on"]("unhandledRejection",function(reason){throw reason});quit_=(status,toThrow)=>{if(keepRuntimeAlive()){process["exitCode"]=status;throw toThrow}logExceptionOnExit(toThrow);process["exit"](status)};Module["inspect"]=function(){return"[Emscripten Module object]"}}else if(ENVIRONMENT_IS_WEB||ENVIRONMENT_IS_WORKER){if(ENVIRONMENT_IS_WORKER){scriptDirectory=self.location.href}else if(typeof document!="undefined"&&document.currentScript){scriptDirectory=document.currentScript.src}if(_scriptDir){scriptDirectory=_scriptDir}if(scriptDirectory.indexOf("blob:")!==0){scriptDirectory=scriptDirectory.substr(0,scriptDirectory.replace(/[?#].*/,"").lastIndexOf("/")+1)}else{scriptDirectory=""}{read_=url=>{var xhr=new XMLHttpRequest;xhr.open("GET",url,false);xhr.send(null);return xhr.responseText};if(ENVIRONMENT_IS_WORKER){readBinary=url=>{var xhr=new XMLHttpRequest;xhr.open("GET",url,false);xhr.responseType="arraybuffer";xhr.send(null);return new Uint8Array(xhr.response)}}readAsync=(url,onload,onerror)=>{var xhr=new XMLHttpRequest;xhr.open("GET",url,true);xhr.responseType="arraybuffer";xhr.onload=()=>{if(xhr.status==200||xhr.status==0&&xhr.response){onload(xhr.response);return}onerror()};xhr.onerror=onerror;xhr.send(null)}}setWindowTitle=title=>document.title=title}else{}var out=Module["print"]||console.log.bind(console);var err=Module["printErr"]||console.warn.bind(console);Object.assign(Module,moduleOverrides);moduleOverrides=null;if(Module["arguments"])arguments_=Module["arguments"];if(Module["thisProgram"])thisProgram=Module["thisProgram"];if(Module["quit"])quit_=Module["quit"];var wasmBinary;if(Module["wasmBinary"])wasmBinary=Module["wasmBinary"];var noExitRuntime=Module["noExitRuntime"]||true;if(typeof WebAssembly!="object"){abort("no native wasm support detected")}var wasmMemory;var ABORT=false;var EXITSTATUS;var buffer,HEAP8,HEAPU8,HEAP16,HEAPU16,HEAP32,HEAPU32,HEAPF32,HEAPF64;function updateGlobalBufferAndViews(buf){buffer=buf;Module["HEAP8"]=HEAP8=new Int8Array(buf);Module["HEAP16"]=HEAP16=new Int16Array(buf);Module["HEAP32"]=HEAP32=new Int32Array(buf);Module["HEAPU8"]=HEAPU8=new Uint8Array(buf);Module["HEAPU16"]=HEAPU16=new Uint16Array(buf);Module["HEAPU32"]=HEAPU32=new Uint32Array(buf);Module["HEAPF32"]=HEAPF32=new Float32Array(buf);Module["HEAPF64"]=HEAPF64=new Float64Array(buf)}var INITIAL_MEMORY=Module["INITIAL_MEMORY"]||20971520;var wasmTable;var __ATPRERUN__=[];var __ATINIT__=[];var __ATPOSTRUN__=[];var runtimeInitialized=false;function keepRuntimeAlive(){return noExitRuntime}function preRun(){if(Module["preRun"]){if(typeof Module["preRun"]=="function")Module["preRun"]=[Module["preRun"]];while(Module["preRun"].length){addOnPreRun(Module["preRun"].shift())}}callRuntimeCallbacks(__ATPRERUN__)}function initRuntime(){runtimeInitialized=true;callRuntimeCallbacks(__ATINIT__)}function postRun(){if(Module["postRun"]){if(typeof Module["postRun"]=="function")Module["postRun"]=[Module["postRun"]];while(Module["postRun"].length){addOnPostRun(Module["postRun"].shift())}}callRuntimeCallbacks(__ATPOSTRUN__)}function addOnPreRun(cb){__ATPRERUN__.unshift(cb)}function addOnInit(cb){__ATINIT__.unshift(cb)}function addOnPostRun(cb){__ATPOSTRUN__.unshift(cb)}var runDependencies=0;var runDependencyWatcher=null;var dependenciesFulfilled=null;function addRunDependency(id){runDependencies++;if(Module["monitorRunDependencies"]){Module["monitorRunDependencies"](runDependencies)}}function removeRunDependency(id){runDependencies--;if(Module["monitorRunDependencies"]){Module["monitorRunDependencies"](runDependencies)}if(runDependencies==0){if(runDependencyWatcher!==null){clearInterval(runDependencyWatcher);runDependencyWatcher=null}if(dependenciesFulfilled){var callback=dependenciesFulfilled;dependenciesFulfilled=null;callback()}}}function abort(what){{if(Module["onAbort"]){Module["onAbort"](what)}}what="Aborted("+what+")";err(what);ABORT=true;EXITSTATUS=1;what+=". Build with -sASSERTIONS for more info.";var e=new WebAssembly.RuntimeError(what);readyPromiseReject(e);throw e}var dataURIPrefix="data:application/octet-stream;base64,";function isDataURI(filename){return filename.startsWith(dataURIPrefix)}function isFileURI(filename){return filename.startsWith("file://")}var wasmBinaryFile;if(Module["locateFile"]){wasmBinaryFile="lib_wasm.wasm";if(!isDataURI(wasmBinaryFile)){wasmBinaryFile=locateFile(wasmBinaryFile)}}else{wasmBinaryFile=new URL("lib_wasm.wasm",import.meta.url).toString()}function getBinary(file){try{if(file==wasmBinaryFile&&wasmBinary){return new Uint8Array(wasmBinary)}if(readBinary){return readBinary(file)}else{throw"both async and sync fetching of the wasm failed"}}catch(err){abort(err)}}function getBinaryPromise(){if(!wasmBinary&&(ENVIRONMENT_IS_WEB||ENVIRONMENT_IS_WORKER)){if(typeof fetch=="function"&&!isFileURI(wasmBinaryFile)){return fetch(wasmBinaryFile,{credentials:"same-origin"}).then(function(response){if(!response["ok"]){throw"failed to load wasm binary file at '"+wasmBinaryFile+"'"}return response["arrayBuffer"]()}).catch(function(){return getBinary(wasmBinaryFile)})}else{if(readAsync){return new Promise(function(resolve,reject){readAsync(wasmBinaryFile,function(response){resolve(new Uint8Array(response))},reject)})}}}return Promise.resolve().then(function(){return getBinary(wasmBinaryFile)})}function createWasm(){var info={"a":asmLibraryArg};function receiveInstance(instance,module){var exports=instance.exports;Module["asm"]=exports;wasmMemory=Module["asm"]["b"];updateGlobalBufferAndViews(wasmMemory.buffer);wasmTable=Module["asm"]["F"];addOnInit(Module["asm"]["c"]);removeRunDependency("wasm-instantiate")}addRunDependency("wasm-instantiate");function receiveInstantiationResult(result){receiveInstance(result["instance"])}function instantiateArrayBuffer(receiver){return getBinaryPromise().then(function(binary){return WebAssembly.instantiate(binary,info)}).then(function(instance){return instance}).then(receiver,function(reason){err("failed to asynchronously prepare wasm: "+reason);abort(reason)})}function instantiateAsync(){if(!wasmBinary&&typeof WebAssembly.instantiateStreaming=="function"&&!isDataURI(wasmBinaryFile)&&!isFileURI(wasmBinaryFile)&&!ENVIRONMENT_IS_NODE&&typeof fetch=="function"){return fetch(wasmBinaryFile,{credentials:"same-origin"}).then(function(response){var result=WebAssembly.instantiateStreaming(response,info);return result.then(receiveInstantiationResult,function(reason){err("wasm streaming compile failed: "+reason);err("falling back to ArrayBuffer instantiation");return instantiateArrayBuffer(receiveInstantiationResult)})})}else{return instantiateArrayBuffer(receiveInstantiationResult)}}if(Module["instantiateWasm"]){try{var exports=Module["instantiateWasm"](info,receiveInstance);return exports}catch(e){err("Module.instantiateWasm callback failed with error: "+e);return false}}instantiateAsync().catch(readyPromiseReject);return{}}function callRuntimeCallbacks(callbacks){while(callbacks.length>0){var callback=callbacks.shift();if(typeof callback=="function"){callback(Module);continue}var func=callback.func;if(typeof func=="number"){if(callback.arg===undefined){getWasmTableEntry(func)()}else{getWasmTableEntry(func)(callback.arg)}}else{func(callback.arg===undefined?null:callback.arg)}}}var wasmTableMirror=[];function getWasmTableEntry(funcPtr){var func=wasmTableMirror[funcPtr];if(!func){if(funcPtr>=wasmTableMirror.length)wasmTableMirror.length=funcPtr+1;wasmTableMirror[funcPtr]=func=wasmTable.get(funcPtr)}return func}function abortOnCannotGrowMemory(requestedSize){abort("OOM")}function _emscripten_resize_heap(requestedSize){var oldSize=HEAPU8.length;requestedSize=requestedSize>>>0;abortOnCannotGrowMemory(requestedSize)}var asmLibraryArg={"a":_emscripten_resize_heap};var asm=createWasm();var ___wasm_call_ctors=Module["___wasm_call_ctors"]=function(){return(___wasm_call_ctors=Module["___wasm_call_ctors"]=Module["asm"]["c"]).apply(null,arguments)};var _alloc_float32array=Module["_alloc_float32array"]=function(){return(_alloc_float32array=Module["_alloc_float32array"]=Module["asm"]["d"]).apply(null,arguments)};var _alloc_uint8array=Module["_alloc_uint8array"]=function(){return(_alloc_uint8array=Module["_alloc_uint8array"]=Module["asm"]["e"]).apply(null,arguments)};var _alloc_uint16array=Module["_alloc_uint16array"]=function(){return(_alloc_uint16array=Module["_alloc_uint16array"]=Module["asm"]["f"]).apply(null,arguments)};var _alloc_uint32array=Module["_alloc_uint32array"]=function(){return(_alloc_uint32array=Module["_alloc_uint32array"]=Module["asm"]["g"]).apply(null,arguments)};var _fill_arr_pt_cloud_3d_sphere=Module["_fill_arr_pt_cloud_3d_sphere"]=function(){return(_fill_arr_pt_cloud_3d_sphere=Module["_fill_arr_pt_cloud_3d_sphere"]=Module["asm"]["h"]).apply(null,arguments)};var _fill_arr_pt_cloud_2d_cube=Module["_fill_arr_pt_cloud_2d_cube"]=function(){return(_fill_arr_pt_cloud_2d_cube=Module["_fill_arr_pt_cloud_2d_cube"]=Module["asm"]["i"]).apply(null,arguments)};var _fill_arr_pt_cloud_3d_cube=Module["_fill_arr_pt_cloud_3d_cube"]=function(){return(_fill_arr_pt_cloud_3d_cube=Module["_fill_arr_pt_cloud_3d_cube"]=Module["asm"]["j"]).apply(null,arguments)};var _move_pts_2d=Module["_move_pts_2d"]=function(){return(_move_pts_2d=Module["_move_pts_2d"]=Module["asm"]["k"]).apply(null,arguments)};var _move_pts_3d=Module["_move_pts_3d"]=function(){return(_move_pts_3d=Module["_move_pts_3d"]=Module["asm"]["l"]).apply(null,arguments)};var _modify_vels_via_curl_noise_2d=Module["_modify_vels_via_curl_noise_2d"]=function(){return(_modify_vels_via_curl_noise_2d=Module["_modify_vels_via_curl_noise_2d"]=Module["asm"]["m"]).apply(null,arguments)};var _compute_hash=Module["_compute_hash"]=function(){return(_compute_hash=Module["_compute_hash"]=Module["asm"]["n"]).apply(null,arguments)};var _draw_pts_2d=Module["_draw_pts_2d"]=function(){return(_draw_pts_2d=Module["_draw_pts_2d"]=Module["asm"]["o"]).apply(null,arguments)};var _draw_surface_pts_2d=Module["_draw_surface_pts_2d"]=function(){return(_draw_surface_pts_2d=Module["_draw_surface_pts_2d"]=Module["asm"]["p"]).apply(null,arguments)};var _draw_pts_3d=Module["_draw_pts_3d"]=function(){return(_draw_pts_3d=Module["_draw_pts_3d"]=Module["asm"]["q"]).apply(null,arguments)};var _pbd_pre=Module["_pbd_pre"]=function(){return(_pbd_pre=Module["_pbd_pre"]=Module["asm"]["r"]).apply(null,arguments)};var _pbd_post=Module["_pbd_post"]=function(){return(_pbd_post=Module["_pbd_post"]=Module["asm"]["s"]).apply(null,arguments)};var _solve_pairwise_constraints_2d=Module["_solve_pairwise_constraints_2d"]=function(){return(_solve_pairwise_constraints_2d=Module["_solve_pairwise_constraints_2d"]=Module["asm"]["t"]).apply(null,arguments)};var _solve_pairwise_constraints_3d=Module["_solve_pairwise_constraints_3d"]=function(){return(_solve_pairwise_constraints_3d=Module["_solve_pairwise_constraints_3d"]=Module["asm"]["u"]).apply(null,arguments)};var _solve_neohookean_constraints_tris_2d=Module["_solve_neohookean_constraints_tris_2d"]=function(){return(_solve_neohookean_constraints_tris_2d=Module["_solve_neohookean_constraints_tris_2d"]=Module["asm"]["v"]).apply(null,arguments)};var _solve_length_constraints_tris_2d=Module["_solve_length_constraints_tris_2d"]=function(){return(_solve_length_constraints_tris_2d=Module["_solve_length_constraints_tris_2d"]=Module["asm"]["w"]).apply(null,arguments)};var _solve_length_constraints_tris_3d=Module["_solve_length_constraints_tris_3d"]=function(){return(_solve_length_constraints_tris_3d=Module["_solve_length_constraints_tris_3d"]=Module["asm"]["x"]).apply(null,arguments)};var _solve_area_constraints_tris_2d=Module["_solve_area_constraints_tris_2d"]=function(){return(_solve_area_constraints_tris_2d=Module["_solve_area_constraints_tris_2d"]=Module["asm"]["y"]).apply(null,arguments)};var _solve_surface_area_constraints_tris_3d=Module["_solve_surface_area_constraints_tris_3d"]=function(){return(_solve_surface_area_constraints_tris_3d=Module["_solve_surface_area_constraints_tris_3d"]=Module["asm"]["z"]).apply(null,arguments)};var _solve_volume_constraints_tris_3d=Module["_solve_volume_constraints_tris_3d"]=function(){return(_solve_volume_constraints_tris_3d=Module["_solve_volume_constraints_tris_3d"]=Module["asm"]["A"]).apply(null,arguments)};var _apply_energy_based_breakage_tris_2d=Module["_apply_energy_based_breakage_tris_2d"]=function(){return(_apply_energy_based_breakage_tris_2d=Module["_apply_energy_based_breakage_tris_2d"]=Module["asm"]["B"]).apply(null,arguments)};var _apply_distance_based_breakage_tris_2d=Module["_apply_distance_based_breakage_tris_2d"]=function(){return(_apply_distance_based_breakage_tris_2d=Module["_apply_distance_based_breakage_tris_2d"]=Module["asm"]["C"]).apply(null,arguments)};var _apply_distance_based_breakage_tris_3d=Module["_apply_distance_based_breakage_tris_3d"]=function(){return(_apply_distance_based_breakage_tris_3d=Module["_apply_distance_based_breakage_tris_3d"]=Module["asm"]["D"]).apply(null,arguments)};var _test=Module["_test"]=function(){return(_test=Module["_test"]=Module["asm"]["E"]).apply(null,arguments)};var calledRun;function ExitStatus(status){this.name="ExitStatus";this.message="Program terminated with exit("+status+")";this.status=status}dependenciesFulfilled=function runCaller(){if(!calledRun)run();if(!calledRun)dependenciesFulfilled=runCaller};function run(args){args=args||arguments_;if(runDependencies>0){return}preRun();if(runDependencies>0){return}function doRun(){if(calledRun)return;calledRun=true;Module["calledRun"]=true;if(ABORT)return;initRuntime();readyPromiseResolve(Module);if(Module["onRuntimeInitialized"])Module["onRuntimeInitialized"]();postRun()}if(Module["setStatus"]){Module["setStatus"]("Running...");setTimeout(function(){setTimeout(function(){Module["setStatus"]("")},1);doRun()},1)}else{doRun()}}Module["run"]=run;if(Module["preInit"]){if(typeof Module["preInit"]=="function")Module["preInit"]=[Module["preInit"]];while(Module["preInit"].length>0){Module["preInit"].pop()()}}run();


  return wasm_.ready
}
);
})();
export default wasm_;