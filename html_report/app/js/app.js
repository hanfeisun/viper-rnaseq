'use strict';

/* App Module */

var viperApp = angular.module('viperApp', [
  'ui.router',
  'viperAppControllers'
]);


viperApp.config(['$stateProvider', '$urlRouterProvider', '$locationProvider', 
        function($stateProvider, $urlRouterProvider, $locationProvider) {

    $locationProvider.html5Mode({
        enabled: true,
        requireBase: false
    });    

    $urlRouterProvider.otherwise('/');

    $stateProvider

        .state('/', {
            url: '/',
            templateUrl: '/app/partials/full-report.html'
        })
        .state('/home', {
            url: '/home',
            templateUrl: '/app/partials/test.html'
        })
        .state('/about',{
            url: '/about',
            templateUrl: '/app/partials/test.html'
        })
        .state('/contact',{
            url: '/contact',
            templateUrl: '/app/partials/test.html'
        })
        .state('/Alignment-Report', {
            url: '/Alignment-Report',
            templateUrl: '/app/partials/section-report.html'
        });
}]);


