'use strict';

/* App Module */

var rnaApp = angular.module('rnaApp', [
  'ngRoute',
  'rnaAppControllers'
]);

rnaApp.config(['$routeProvider',
  function($routeProvider) {
    $routeProvider.
      when('/full-report', {
        templateUrl: 'partials/full-report.html',
        controller: 'FullReportCtrl'
      }).
      when('/full-report/:sectionId', {
        templateUrl: 'partials/section-report.html',
        controller: 'SectionReportCtrl'
      }).
      otherwise({
        redirectTo: '/full-report'
      });
  }]);
