var gulp = require('gulp');

var concat = require('gulp-concat');
var stripDebug = require('gulp-strip-debug');
var uglify = require('gulp-uglify');
var changed = require('gulp-changed');
var merge = require('gulp-merge');
var pump = require('pump');
  // var gzip = require('gulp-gzip');

// js
gulp.task('scripts.min', function(cb){
  pump([
      gulp.src([
         'wrap/on.js',
         'libs/numeric.core.js',
         'libs/jacobi.js',
         'wrap/off.js'
      ]),
      concat('numerico.min.js'),
      stripDebug(),
      uglify(),
      gulp.dest('.')
    ],
    cb
  );
});


gulp.task('scripts', function(cb){
  pump([
      gulp.src([
         'wrap/on.js',
         'libs/numeric.core.js',
         'libs/jacobi.js',
         'wrap/off.js'
      ]),
      concat('numerico.js'),
      gulp.dest('.')
    ],
    cb
  );
});


// watch
gulp.task('watch', function(){
   gulp.watch(['libs/*.js'], ['scripts', 'scripts.min']);
});


// default gulp task
gulp.task('default', ['scripts', 'scripts.min', 'watch']);

